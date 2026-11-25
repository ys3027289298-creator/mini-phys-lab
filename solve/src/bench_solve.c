/*
 * bench_solve: 30-sphere sequential-impulse contact solver + profiler.
 *
 * This file reconstructs the "30 balls, ~8 ms contact solve" workload described
 * in the task, because the repository only shipped single-contact self tests.
 *
 * Timed contact solve = contact enumeration + per-contact derived-data setup +
 * Gauss-Seidel impulse iterations.  Two code paths are provided:
 *
 *   legacy : fixed iteration count K, derived data (normal unitization,
 *            tangents, effective masses, Baumgarte bias) recomputed for every
 *            contact on every iteration -- the "repeated per-iteration
 *            computation" pathology, with separate timers for recompute vs the
 *            impulse apply core.
 *   cached : derived data computed once after enumeration + convergence stop
 *            when the max normal accumulated-impulse change between adjacent
 *            iterations drops below epsilon.
 *
 * Modes:
 *   bench_solve                         final verification (default eps)
 *   bench_solve baseline [K] [reps]     three-way phase timing
 *   bench_solve itersweep [reps]        CSV K vs time (legacy + cached)
 *   bench_solve converge [Kmax]         CSV per-iteration impulse delta
 *   bench_solve sim <eps> [frames]      one JSON line: residual + timing
 *   bench_solve final <eps> [reps] [frames]
 *   bench_solve scenario                contact-count / geometry summary
 *
 * All timing uses clock_gettime(CLOCK_MONOTONIC); nothing is estimated.
 */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* ------------------------------------------------------------------ V3 -- */

typedef struct { float x, y, z; } V3;

static inline V3 v3(float x, float y, float z) { V3 r = {x, y, z}; return r; }
static inline V3 vadd(V3 a, V3 b) { return v3(a.x+b.x, a.y+b.y, a.z+b.z); }
static inline V3 vsub(V3 a, V3 b) { return v3(a.x-b.x, a.y-b.y, a.z-b.z); }
static inline V3 vmul(V3 a, float s) { return v3(a.x*s, a.y*s, a.z*s); }
static inline float vdot(V3 a, V3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline float vlen(V3 a) { return sqrtf(vdot(a, a)); }

/* --------------------------------------------------------------- world --- */

#define NBALLS   30
#define NPLANES  5
#define NBODIES  (NBALLS + NPLANES)
#define MAXC     768
#define MAXIT    8192

/* chosen by artifacts/eps_sweep.csv; see artifacts/contact_solve_report.md.
 * 0.03 is the loosest deterministic pass (res 9.87e-4, only 1.3% margin);
 * 0.01 buys 4.5x residual margin at one extra iteration (~6 us). */
#define DEFAULT_EPS 0.01f
/* legacy fixed iteration count. 4096 is the historical "crank it until stable"
 * setting; on the reference Ryzen box its uninstrumented median is ~8 ms for
 * the 30-ball compacted frame, matching the reported profile. */
#define LEGACY_ITERS 4096

typedef struct {
    V3 p, v;
    float r;
    float im;                 /* 1/mass; 0 for static plane proxies */
} Body;

typedef struct {
    V3 p;                     /* a point on the plane */
    V3 n;                     /* inward unit normal */
} Plane;

typedef struct {
    int a, b;                 /* b indexes a static proxy body (im == 0) */
    V3 n;                     /* from b to a; +impulse pushes a along n */
    float depth;              /* positive = penetration */

    /* ---- derived invariant data (the repeated-computation item) ---- */
    V3 t1, t2;
    float kn;                 /* normal effective mass 1/(im_a+im_b) */
    float kt;                 /* tangent effective mass (same for spheres) */
    float bias;               /* Baumgarte target normal speed */

    /* ---- accumulated impulses ---- */
    float pn, pt1, pt2;
} Contact;

typedef struct {
    Body bd[NBODIES];
    Plane pl[NPLANES];
    Contact ct[MAXC];
    int nc;

    float dt;
    float gravity;
    float beta;               /* Baumgarte position correction fraction */
    float mu;                 /* Coulomb friction coefficient */
    float skin;               /* speculative contact skin */
} World;

/* ------------------------------------------------------------- timing ---- */

static uint64_t now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ull + (uint64_t)ts.tv_nsec;
}

static int cmp_u64(const void* A, const void* B) {
    uint64_t a = *(const uint64_t*)A, b = *(const uint64_t*)B;
    return (a > b) - (a < b);
}

static int cmp_int(const void* A, const void* B) {
    int a = *(const int*)A, b = *(const int*)B;
    return (a > b) - (a < b);
}

typedef struct {
    uint64_t med, p95, mn, mx;
} Dist;

static Dist summarize(uint64_t* x, int n) {
    qsort(x, (size_t)n, sizeof(x[0]), cmp_u64);
    Dist d;
    d.med = x[n / 2];
    d.p95 = x[(int)((n - 1) * 0.95)];
    d.mn = x[0];
    d.mx = x[n - 1];
    return d;
}

/* prevent the optimizer from deleting the measured work */
static float g_sink = 0.0f;

/* ------------------------------------------------------------ scenario --- */

static unsigned hsh(int i, int salt) {
    unsigned x = (unsigned)(i * 73856093) ^ (unsigned)(salt * 19349663);
    x ^= x >> 13; x *= 1274126177u; x ^= x >> 16;
    return x;
}
static float u01(int i, int salt) { return (hsh(i, salt) & 0xffffff) / (float)0x1000000; }

/*
 * 5 x 3 x 2 cubic close grid = 30 unit spheres, floor + 4 walls.
 * overlap > 0 shrinks the grid spacing to pre-load every contact (stress
 * frame); overlap == 0 puts neighbours exactly touching (residual runs).
 */
static void init_pile(World* w, float overlap, float down_speed, int perturb) {
    memset(w, 0, sizeof(*w));
    w->dt = 1.0f / 120.0f;
    w->gravity = 9.81f;
    w->beta = 0.2f;
    w->mu = 0.3f;
    w->skin = 0.005f;

    const float R = 1.0f;
    const float S = 2.0f * R - overlap;
    const float WX = 2.0f * S + R;     /* +/- 4.98 when overlap .01 */
    const float WZ = 1.0f * S + R;     /* +/- 2.99 */

    /* static proxy bodies first (indices NBALLS..NBODIES-1) */
    for (int i = NBALLS; i < NBODIES; i++) {
        w->bd[i].im = 0.0f;
        w->bd[i].r = 0.0f;
    }
    w->pl[0] = (Plane){ v3(0, 0, 0),    v3(0, 1, 0) };   /* floor */
    w->pl[1] = (Plane){ v3(-WX, 0, 0),  v3(1, 0, 0) };
    w->pl[2] = (Plane){ v3( WX, 0, 0),  v3(-1, 0, 0) };
    w->pl[3] = (Plane){ v3(0, 0, -WZ),  v3(0, 0, 1) };
    w->pl[4] = (Plane){ v3(0, 0,  WZ),  v3(0, 0, -1) };

    int k = 0;
    for (int ly = 0; ly < 2; ly++)
    for (int iz = 0; iz < 3; iz++)
    for (int ix = 0; ix < 5; ix++) {
        Body* b = &w->bd[k];
        b->r = R;
        b->im = 1.0f;                 /* mass 1 */
        b->p = v3((ix - 2) * S, R + ly * S, (iz - 1) * S);
        b->v = v3(0.0f, -down_speed, 0.0f);
        if (perturb) {
            /* velocity-only symmetry breaking: positional jitter would inject
             * up to ~1e-3 initial penetration and confound the residual curve */
            b->v.x = (u01(k, 1) - 0.5f) * 0.20f;
            b->v.z = (u01(k, 2) - 0.5f) * 0.20f;
        }
        k++;
    }
}

/* ---------------------------------------------------------- enumeration -- */

static void add_contact(World* w, int a, int b, V3 n, float depth) {
    if (w->nc >= MAXC) { fprintf(stderr, "contact overflow\n"); exit(2); }
    Contact* c = &w->ct[w->nc++];
    memset(c, 0, sizeof(*c));
    c->a = a; c->b = b; c->n = n; c->depth = depth;
}

static void enumerate(World* w) {
    w->nc = 0;
    const float skin = w->skin;

    for (int i = 0; i < NBALLS; i++) {
        Body* A = &w->bd[i];
        for (int j = i + 1; j < NBALLS; j++) {
            Body* B = &w->bd[j];
            V3 d = vsub(A->p, B->p);
            float rs = A->r + B->r;
            float lim = rs + skin;
            float d2 = vdot(d, d);
            if (d2 > lim * lim) continue;
            float dist = sqrtf(d2);
            V3 n;
            if (dist > 1e-9f) n = vmul(d, 1.0f / dist);
            else n = v3(0, 1, 0);
            add_contact(w, i, j, n, rs - dist);
        }
        for (int q = 0; q < NPLANES; q++) {
            V3 rel = vsub(A->p, w->pl[q].p);
            float dist = vdot(rel, w->pl[q].n);   /* signed gap to plane */
            if (dist > A->r + skin) continue;
            add_contact(w, i, NBALLS + q, w->pl[q].n, A->r - dist);
        }
    }
}

/* ----------------------------------------- derived data (recomputed) ----- */

/*
 * Everything here is invariant across velocity iterations for fixed positions:
 * normal unitization, orthonormal friction basis, effective masses and the
 * Baumgarte bias.  A legacy solver recomputes it for every contact on every
 * iteration; the cached path calls it once.
 */
static void contact_prepare(World* w, Contact* c) {
    float l = vlen(c->n);
    if (l > 0.0f) c->n = vmul(c->n, 1.0f / l);

    V3 axis = (fabsf(c->n.y) < 0.95f) ? v3(0, 1, 0) : v3(1, 0, 0);
    V3 t1 = v3(
        axis.y * c->n.z - axis.z * c->n.y,
        axis.z * c->n.x - axis.x * c->n.z,
        axis.x * c->n.y - axis.y * c->n.x);
    float tl = vlen(t1);
    if (tl > 1e-9f) t1 = vmul(t1, 1.0f / tl);
    else t1 = v3(1, 0, 0);
    V3 t2 = v3(
        c->n.y * t1.z - c->n.z * t1.y,
        c->n.z * t1.x - c->n.x * t1.z,
        c->n.x * t1.y - c->n.y * t1.x);
    c->t1 = t1; c->t2 = t2;

    float wsum = w->bd[c->a].im + w->bd[c->b].im;
    c->kn = (wsum > 0.0f) ? 1.0f / wsum : 0.0f;
    c->kt = c->kn;   /* spherical particles: same inverse mass along tangent */

    float pen = c->depth;
    c->bias = (pen > 0.0f) ? (w->beta / w->dt) * pen : 0.0f;
}

static void prepare_all(World* w) {
    for (int i = 0; i < w->nc; i++) contact_prepare(w, &w->ct[i]);
}

/* ------------------------------------------------- impulse application --- */

static inline void apply_vec_impulse(World* w, const Contact* c, V3 j) {
    Body* a = &w->bd[c->a];
    Body* b = &w->bd[c->b];
    a->v = vadd(a->v, vmul(j, a->im));
    b->v = vsub(b->v, vmul(j, b->im));   /* static proxy has im == 0 */
}

/* returns |delta accumulated normal impulse| */
static inline float solve_normal(World* w, Contact* c) {
    Body* a = &w->bd[c->a];
    Body* b = &w->bd[c->b];
    float vn = vdot(vsub(a->v, b->v), c->n);
    float dpn = c->kn * (c->bias - vn);
    float old = c->pn;
    float pn = old + dpn;
    if (pn < 0.0f) pn = 0.0f;
    float d = pn - old;
    c->pn = pn;
    if (d != 0.0f) apply_vec_impulse(w, c, vmul(c->n, d));
    return fabsf(d);
}

/* boxed Coulomb friction along two tangent axes; returns max |d tangent P| */
static inline float solve_friction(World* w, Contact* c) {
    Body* a = &w->bd[c->a];
    Body* b = &w->bd[c->b];
    float lim = w->mu * c->pn;
    float dmax = 0.0f;

    V3 rv = vsub(a->v, b->v);
    float vt1 = vdot(rv, c->t1);
    float old1 = c->pt1;
    float n1 = old1 - vt1 * c->kt;
    if (n1 >  lim) n1 =  lim;
    if (n1 < -lim) n1 = -lim;
    float d1 = n1 - old1;
    c->pt1 = n1;
    if (d1 != 0.0f) apply_vec_impulse(w, c, vmul(c->t1, d1));
    if (fabsf(d1) > dmax) dmax = fabsf(d1);

    rv = vsub(a->v, b->v);
    float vt2 = vdot(rv, c->t2);
    float old2 = c->pt2;
    float n2 = old2 - vt2 * c->kt;
    if (n2 >  lim) n2 =  lim;
    if (n2 < -lim) n2 = -lim;
    float d2 = n2 - old2;
    c->pt2 = n2;
    if (d2 != 0.0f) apply_vec_impulse(w, c, vmul(c->t2, d2));
    if (fabsf(d2) > dmax) dmax = fabsf(d2);
    return dmax;
}

/* ------------------------------------------- solver variants / timing --- */

typedef struct {
    uint64_t en, redund, core;   /* accumulated ns across a batch */
    int iters;
} PhaseTime;

/* legacy: fixed K; every iteration rebuilds all derived data, then applies */
static void solve_legacy(World* w, int K, int instrument, PhaseTime* pt) {
    uint64_t t0 = now_ns();
    enumerate(w);
    uint64_t t1 = now_ns();
    if (instrument && pt) pt->en += t1 - t0;

    for (int k = 0; k < K; k++) {
        if (instrument) t0 = now_ns();
        for (int i = 0; i < w->nc; i++) contact_prepare(w, &w->ct[i]);
        if (instrument) { t1 = now_ns(); pt->redund += t1 - t0; t0 = t1; }

        for (int i = 0; i < w->nc; i++) {
            float d = solve_normal(w, &w->ct[i]);
            solve_friction(w, &w->ct[i]);
            g_sink += d * 1e-12f;
        }
        if (instrument) { t1 = now_ns(); pt->core += t1 - t0; }
    }
    if (pt) pt->iters = K;
}

/* cached derived data + epsilon convergence stop on normal impulse delta */
static int solve_cached(World* w, float eps, int maxit, float* dlog) {
    enumerate(w);
    prepare_all(w);
    int used = 0;
    for (int k = 0; k < maxit; k++) {
        float dmax = 0.0f;
        for (int i = 0; i < w->nc; i++) {
            float dn = solve_normal(w, &w->ct[i]);
            float df = solve_friction(w, &w->ct[i]);
            if (dn > dmax) dmax = dn;
            g_sink += (dn + df) * 1e-12f;
        }
        used = k + 1;
        if (dlog) dlog[k] = dmax;
        if (dmax < eps) break;
    }
    return used;
}

/* position update only. Gravity/forces are applied exactly once per frame
 * BEFORE the solve (semi-implicit Euler: v += g*dt; solve(v); x += v*dt).
 * Applying gravity here too would inject an unsolved post-solve velocity and
 * create a spurious steady penetration g*dt^2/beta. */
static void integrate(World* w) {
    for (int i = 0; i < NBALLS; i++) {
        Body* b = &w->bd[i];
        b->p = vadd(b->p, vmul(b->v, w->dt));
    }
}

/* independent O(N^2) ground-truth penetration check, no skin */
static float max_penetration(const World* w) {
    float m = 0.0f;
    for (int i = 0; i < NBALLS; i++) {
        const Body* a = &w->bd[i];
        for (int j = i + 1; j < NBALLS; j++) {
            const Body* b = &w->bd[j];
            V3 d = vsub(a->p, b->p);
            float dist = sqrtf(vdot(d, d));
            float pen = a->r + b->r - dist;
            if (pen > m) m = pen;
        }
        for (int q = 0; q < NPLANES; q++) {
            float gap = vdot(vsub(a->p, w->pl[q].p), w->pl[q].n);
            float pen = a->r - gap;
            if (pen > m) m = pen;
        }
    }
    return m;
}

/* ------------------------------------------------------------- modes ----- */

/* cold stress-frame snapshot: a gravity-compacted, feasible dense contact web
 * in its post-gravity / pre-solve state, with zero accumulated impulse. */
static World g_stress;

static void make_stress(void) {
    init_pile(&g_stress, 0.0f, 0.0f, 1);
    /* compact under gravity with a tight reference solve */
    for (int f = 0; f < 24; f++) {
        for (int i = 0; i < NBALLS; i++) g_stress.bd[i].v.y -= g_stress.gravity * g_stress.dt;
        solve_cached(&g_stress, 1e-6f, 4096, NULL);
        integrate(&g_stress);
    }
    /* next frame's gravity load is already present, as in run_sim pre-solve */
    for (int i = 0; i < NBALLS; i++) g_stress.bd[i].v.y -= g_stress.gravity * g_stress.dt;
}

static int mode_scenario(void) {
    World w;
    init_pile(&w, 0.01f, 3.0f, 0);
    enumerate(&w);
    int bb = 0, pl = 0;
    for (int i = 0; i < w.nc; i++) if (w.ct[i].b < NBALLS) bb++; else pl++;
    printf("balls=%d bodies=%d contacts=%d ball_ball=%d ball_plane=%d "
           "max_pairs=%d dt=%.6f beta=%.2f mu=%.2f\n",
           NBALLS, NBODIES, w.nc, bb, pl,
           NBALLS * (NBALLS - 1) / 2 + NBALLS * NPLANES,
           w.dt, w.beta, w.mu);
    World z; init_pile(&z, 0.0f, 0.0f, 1); enumerate(&z);
    printf("touching_scenario_contacts=%d initial_max_pen=%.3e\n",
           z.nc, max_penetration(&z));
    return 0;
}

static int mode_baseline(int K, int reps) {
    if (K <= 0) K = LEGACY_ITERS;
    if (reps <= 0) reps = 500;

    /* pass 1: uninstrumented outer timing, robust medians */
    uint64_t* outer = calloc((size_t)reps, sizeof(uint64_t));
    for (int r = 0; r < reps; r++) {
        World w = g_stress;
        uint64_t t0 = now_ns();
        solve_legacy(&w, K, 0, NULL);
        outer[r] = now_ns() - t0;
        g_sink += w.ct[0].pn * 1e-15f;
    }
    Dist do_ = summarize(outer, reps);

    /* pass 2: per-rep instrumented phase times, also summarised as medians so
     * the three items share one robust statistic (means get dragged by scheduler
     * stalls on a shared host). */
    uint64_t* ae = calloc((size_t)reps, sizeof(uint64_t));
    uint64_t* ar = calloc((size_t)reps, sizeof(uint64_t));
    uint64_t* ac = calloc((size_t)reps, sizeof(uint64_t));
    uint64_t* at = calloc((size_t)reps, sizeof(uint64_t));
    for (int r = 0; r < reps; r++) {
        World w = g_stress;
        PhaseTime pt = {0};
        uint64_t t0 = now_ns();
        solve_legacy(&w, K, 1, &pt);
        at[r] = now_ns() - t0;
        ae[r] = pt.en; ar[r] = pt.redund; ac[r] = pt.core;
    }
    Dist de = summarize(ae, reps), dr = summarize(ar, reps),
         dc = summarize(ac, reps), dt = summarize(at, reps);
    double ren = de.med / 1e3, rrd = dr.med / 1e3, rco = dc.med / 1e3;
    double sum = ren + rrd + rco;

    printf("# legacy fixed K=%d reps=%d contacts=%d (medians)\n", K, reps, g_stress.nc);
    printf("baseline_outer_median_us=%.3f p95_us=%.3f min_us=%.3f max_us=%.3f\n",
           do_.med / 1e3, do_.p95 / 1e3, do_.mn / 1e3, do_.mx / 1e3);
    printf("phase_enum_us=%.3f share=%.4f\n", ren, ren / sum);
    printf("phase_redund_us=%.3f share=%.4f\n", rrd, rrd / sum);
    printf("phase_core_us=%.3f share=%.4f\n", rco, rco / sum);
    printf("phase_sum_us=%.3f instrumented_outer_median_us=%.3f "
           "timer_overhead_vs_uninstrumented=%.2f%%\n",
           sum, dt.med / 1e3,
           100.0 * ((double)dt.med - (double)do_.med) / (double)do_.med);
    printf("JSON {\"mode\":\"baseline\",\"K\":%d,\"reps\":%d,\"contacts\":%d,"
           "\"outer_median_us\":%.3f,\"outer_p95_us\":%.3f,"
           "\"enum_us\":%.3f,\"redund_us\":%.3f,\"core_us\":%.3f,"
           "\"enum_share\":%.4f,\"redund_share\":%.4f,\"core_share\":%.4f,"
           "\"timer_overhead_pct\":%.2f}\n",
           K, reps, g_stress.nc, do_.med / 1e3, do_.p95 / 1e3,
           ren, rrd, rco, ren / sum, rrd / sum, rco / sum,
           100.0 * ((double)dt.med - (double)do_.med) / (double)do_.med);
    free(outer); free(ae); free(ar); free(ac); free(at);
    return 0;
}

static int mode_itersweep(int reps) {
    if (reps <= 0) reps = 200;
    const int Ks[] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096};
    const int NK = (int)(sizeof(Ks) / sizeof(Ks[0]));
    printf("K,legacy_median_us,legacy_p95_us,cachedfix_median_us,per_iter_legacy_ns,per_iter_cached_ns\n");
    for (int q = 0; q < NK; q++) {
        int K = Ks[q];
        uint64_t* a = calloc((size_t)reps, sizeof(uint64_t));
        uint64_t* c = calloc((size_t)reps, sizeof(uint64_t));
        for (int r = 0; r < reps; r++) {
            World w = g_stress;
            uint64_t t0 = now_ns();
            solve_legacy(&w, K, 0, NULL);
            a[r] = now_ns() - t0;
        }
        for (int r = 0; r < reps; r++) {
            World w = g_stress;
            /* eps=0 + maxit=K == fixed K cached path */
            uint64_t t0 = now_ns();
            solve_cached(&w, 0.0f, K, NULL);
            c[r] = now_ns() - t0;
        }
        Dist da = summarize(a, reps), dc = summarize(c, reps);
        double la = da.med / 1e3, lc = dc.med / 1e3;
        printf("%d,%.3f,%.3f,%.3f,%.2f,%.2f\n", K, la, da.p95 / 1e3, lc,
               da.med / (double)K, dc.med / (double)K);
        free(a); free(c);
    }
    return 0;
}

static int mode_converge(int Kmax) {
    if (Kmax <= 0 || Kmax > MAXIT) Kmax = 4096;
    World w = g_stress;
    float* log = calloc((size_t)Kmax, sizeof(float));
    int used = solve_cached(&w, 0.0f, Kmax, log);
    printf("iter,dpn_max\n");
    for (int k = 0; k < used; k++) printf("%d,%.9g\n", k + 1, log[k]);
    free(log);
    return 0;
}

typedef struct {
    float eps;
    int frames;
    float res_max;       /* worst over all frames */
    float res_steady;    /* worst over second half */
    double res_sum; long res_n;
    uint64_t solve_med, solve_p95;
    int it_med, it_max;
    int hit_cap;
} SimResult;

static SimResult run_sim(float eps, int frames, int verbose_timing) {
    SimResult R;
    memset(&R, 0, sizeof(R));
    R.eps = eps; R.frames = frames;

    World w;
    init_pile(&w, 0.0f, 0.0f, 1);   /* exact-touching, mild velocity jitter */

    uint64_t* ts = calloc((size_t)frames, sizeof(uint64_t));
    int* its = calloc((size_t)frames, sizeof(int));
    uint64_t tall = 0;

    for (int f = 0; f < frames; f++) {
        for (int i = 0; i < NBALLS; i++) w.bd[i].v.y -= w.gravity * w.dt;

        uint64_t t0 = now_ns();
        int it = solve_cached(&w, eps, 4096, NULL);
        uint64_t dt_ns = now_ns() - t0;
        tall += dt_ns;
        ts[f] = dt_ns;
        its[f] = it;
        if (it >= 4096) R.hit_cap = 1;

        integrate(&w);

        float pr = max_penetration(&w);
        if (pr > R.res_max) R.res_max = pr;
        if (f >= frames / 2 && pr > R.res_steady) R.res_steady = pr;
        R.res_sum += pr; R.res_n++;
    }

    Dist dt = summarize(ts, frames);
    qsort(its, (size_t)frames, sizeof(int), cmp_int);
    R.solve_med = dt.med; R.solve_p95 = dt.p95;
    R.it_med = its[frames / 2];
    R.it_max = its[frames - 1];
    (void)verbose_timing; (void)tall;
    free(ts); free(its);
    return R;
}

static int mode_sim(float eps, int frames) {
    if (frames <= 0) frames = 480;
    SimResult R = run_sim(eps, frames, 0);
    printf("JSON {\"mode\":\"sim\",\"eps\":%.9g,\"frames\":%d,"
           "\"res_max\":%.6e,\"res_steady\":%.6e,\"res_mean\":%.6e,"
           "\"solve_med_us\":%.3f,\"solve_p95_us\":%.3f,"
           "\"it_med\":%d,\"it_max\":%d,\"hit_cap\":%d}\n",
           eps, frames, R.res_max, R.res_steady, R.res_sum / R.res_n,
           R.solve_med / 1e3, R.solve_p95 / 1e3,
           R.it_med, R.it_max, R.hit_cap);
    return 0;
}

/* robust cold-stress timing for one eps */
static Dist time_final(float eps, int reps, int* it_out) {
    uint64_t* ts = calloc((size_t)reps, sizeof(uint64_t));
    int itlast = 0;
    for (int r = 0; r < reps; r++) {
        World w = g_stress;
        uint64_t t0 = now_ns();
        itlast = solve_cached(&w, eps, 4096, NULL);
        ts[r] = now_ns() - t0;
    }
    if (it_out) *it_out = itlast;
    Dist d = summarize(ts, reps);
    free(ts);
    return d;
}

static int mode_final(float eps, int reps, int frames) {
    if (!(eps > 0.0f)) eps = DEFAULT_EPS;
    if (reps <= 0) reps = 2001;
    if (frames <= 0) frames = 480;

    int itcold;
    Dist td = time_final(eps, reps, &itcold);
    SimResult R = run_sim(eps, frames, 0);

    /* threshold inverse solve bookkeeping numbers (formula checked in report) */
    float wmax = 2.0f;              /* two unit-mass balls */
    float eps_x = 1e-3f;            /* r = 1, so normalized == absolute */
    float eps_lambda_naive = eps_x / (g_stress.dt * wmax);
    /* steady Baumgarte balance: eps_x ~ w*eps_lambda*dt/beta */
    float eps_lambda_baum = g_stress.beta * eps_x / (g_stress.dt * wmax);

    int pass_time = (td.p95 <= 4000000ull);
    int pass_res = (R.res_max < 1e-3f);

    printf("# final verification eps=%.6g reps=%d frames=%d\n", eps, reps, frames);
    printf("cold_stress_iters=%d\n", itcold);
    printf("solve_median_us=%.3f solve_p95_us=%.3f solve_max_us=%.3f "
           "(budget 4000 us)\n", td.med / 1e3, td.p95 / 1e3, td.mx / 1e3);
    printf("residual_max=%.6e residual_steady=%.6e (limit 1.000e-03)\n",
           R.res_max, R.res_steady);
    printf("sim_solve_med_us=%.3f sim_solve_p95_us=%.3f it_med=%d it_max=%d\n",
           R.solve_med / 1e3, R.solve_p95 / 1e3, R.it_med, R.it_max);
    printf("derived_eps_lambda_naive=%.6f (eps_x/(dt*wmax))\n", eps_lambda_naive);
    printf("derived_eps_lambda_baum=%.6f (beta*eps_x/(dt*wmax), q=1 local)\n",
           eps_lambda_baum);
    printf("VERDICT pass_time=%d pass_residual=%d overall=%d\n",
           pass_time, pass_res, pass_time && pass_res);
    printf("JSON {\"mode\":\"final\",\"eps\":%.9g,\"reps\":%d,\"frames\":%d,"
           "\"cold_iters\":%d,\"median_us\":%.3f,\"p95_us\":%.3f,\"max_us\":%.3f,"
           "\"res_max\":%.6e,\"res_steady\":%.6e,"
           "\"sim_med_us\":%.3f,\"sim_p95_us\":%.3f,"
           "\"it_med\":%d,\"it_max\":%d,\"hit_cap\":%d,"
           "\"pass_time\":%d,\"pass_residual\":%d,\"overall\":%d}\n",
           eps, reps, frames, itcold, td.med / 1e3, td.p95 / 1e3, td.mx / 1e3,
           R.res_max, R.res_steady, R.solve_med / 1e3, R.solve_p95 / 1e3,
           R.it_med, R.it_max, R.hit_cap,
           pass_time, pass_res, pass_time && pass_res);
    return (pass_time && pass_res) ? 0 : 1;
}

static int mode_diag(float overlap, float mu, int K) {
    World w;
    init_pile(&w, overlap, 3.0f, 0);
    w.mu = mu;
    enumerate(&w); prepare_all(&w);
    for (int k = 0; k < K; k++) {
        int worst = -1; float wd = -1;
        for (int i = 0; i < w.nc; i++) {
            Contact* c = &w.ct[i];
            float dn = solve_normal(&w, c);
            solve_friction(&w, c);
            if (dn > wd) { wd = dn; worst = i; }
        }
        if (k < 12 || k == 49 || k == 99 || k == 499 || k == 999 || k == K - 1) {
            Contact* c = &w.ct[worst];
            float vn = vdot(vsub(w.bd[c->a].v, w.bd[c->b].v), c->n);
            printf("k=%4d worst=c%03d a=%2d b=%2d dP=%.5f vn=%.5f pn=%.4f "
                   "depth=%.4f n=(%.2f,%.2f,%.2f)\n",
                   k + 1, worst, c->a, c->b, wd, vn, c->pn, c->depth,
                   c->n.x, c->n.y, c->n.z);
        }
    }
    return 0;
}

int main(int argc, char** argv) {
    make_stress();
    const char* mode = (argc > 1) ? argv[1] : "all";

    if (!strcmp(mode, "help") || !strcmp(mode, "-h") || !strcmp(mode, "--help")) {
        printf("usage: bench_solve [all|scenario|baseline [K reps]|"
               "itersweep [reps]|converge [Kmax]|sim eps [frames]|"
               "final eps [reps frames]|selftest|diag ...]\n");
        return 0;
    }

    if (!strcmp(mode, "scenario")) return mode_scenario();
    if (!strcmp(mode, "diag"))
        return mode_diag(argc > 2 ? (float)atof(argv[2]) : 0.01f,
                         argc > 3 ? (float)atof(argv[3]) : 0.3f,
                         argc > 4 ? atoi(argv[4]) : 200);
    if (!strcmp(mode, "baseline"))
        return mode_baseline(argc > 2 ? atoi(argv[2]) : 0,
                             argc > 3 ? atoi(argv[3]) : 300);
    if (!strcmp(mode, "itersweep"))
        return mode_itersweep(argc > 2 ? atoi(argv[2]) : 200);
    if (!strcmp(mode, "converge"))
        return mode_converge(argc > 2 ? atoi(argv[2]) : 4096);
    if (!strcmp(mode, "sim")) {
        float eps = (argc > 2) ? (float)atof(argv[2]) : DEFAULT_EPS;
        return mode_sim(eps, argc > 3 ? atoi(argv[3]) : 480);
    }
    if (!strcmp(mode, "final")) {
        float eps = (argc > 2) ? (float)atof(argv[2]) : DEFAULT_EPS;
        return mode_final(eps, argc > 3 ? atoi(argv[3]) : 2001,
                          argc > 4 ? atoi(argv[4]) : 480);
    }
    if (!strcmp(mode, "selftest")) {
        /* physics sanity checks */
        int fails = 0;
        World w; init_pile(&w, 0.0f, 0.0f, 0);
        enumerate(&w); prepare_all(&w);
        int nc0 = w.nc;
        int it = solve_cached(&w, 1e-7f, 4096, NULL);
        int finite = 1;
        for (int i = 0; i < NBALLS; i++)
            if (!isfinite(w.bd[i].v.x + w.bd[i].v.y + w.bd[i].v.z)) finite = 0;
        float res = max_penetration(&w);
        printf("selftest contacts=%d iters=%d finite=%d pen_after=%.3e\n",
               nc0, it, finite, res);
        if (!(nc0 > 80)) { printf("FAIL contacts too few\n"); fails++; }
        if (!finite) { printf("FAIL non-finite velocity\n"); fails++; }
        printf("selftest_all_ok=%d\n", fails == 0);
        return fails ? 1 : 0;
    }

    if (strcmp(mode, "all") != 0) {
        fprintf(stderr, "unknown mode: %s (try --help)\n", mode);
        return 2;
    }

    /* default "all": full reproducible verdict */
    int rc1 = mode_scenario();
    printf("\n");
    mode_baseline(LEGACY_ITERS, 300);
    printf("\n");
    int rc2 = mode_final(DEFAULT_EPS, 2001, 600);
    (void)rc1;
    return rc2;
}
