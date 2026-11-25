
/* am_bench: generates the data behind the precision decision.
 *   <outdir>/am_sensitivity.csv  rel angular-momentum error vs arm error [mm]
 *   <outdir>/am_tradeoff.csv     per precision tier: max error + ns/collision
 *   <outdir>/am_micro_timing.csv raw cost of one cross product / 2D marginal
 * Run: ./am_bench [outdir]   (default outdir = .)
 */
#define _POSIX_C_SOURCE 199309L
#include "am.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

static double now_s(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static volatile double g_sink;

static int cmp_double(const void* pa, const void* pb) {
    double x = *(const double*)pa, y = *(const double*)pb;
    return (x > y) - (x < y);
}

/* deterministic PRNG so the curves are reproducible */
static unsigned long long st = 0x9e3779b97f4a7c15ULL;
static double urand(void) {
    unsigned long long x = st;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    st = x;
    return (double)((x >> 11) & ((1ULL << 53) - 1)) / (double)(1ULL << 53);
}
static double frand(double lo, double hi) { return lo + (hi - lo) * urand(); }

typedef struct { am_body A, B; am_contact c; } Case;

static void random_case(Case* q) {
    memset(q, 0, sizeof(*q));
    double nx = frand(-1, 1), ny = frand(-1, 1), nz = frand(-1, 1);
    double nl = sqrt(nx * nx + ny * ny + nz * nz) + 1e-12;
    q->c.n = am_v(nx / nl, ny / nl, nz / nl);
    q->c.e = frand(0.0, 1.0);
    q->c.p = am_v(0, 0, 0);
    am_vec3 rA = am_v(frand(-0.6, 0.6), frand(-0.6, 0.6), frand(-0.6, 0.6));
    am_vec3 rB = am_v(frand(-0.6, 0.6), frand(-0.6, 0.6), frand(-0.6, 0.6));
    q->A.x = am_scale(rA, -1.0);
    q->B.x = am_scale(rB, -1.0);
    q->A.m = frand(0.5, 5.0); q->B.m = frand(0.5, 5.0);
    q->A.I = frand(0.02, 0.6); q->B.I = frand(0.02, 0.6);
    q->A.w = am_v(frand(-6, 6), frand(-6, 6), frand(-6, 6));
    q->B.w = am_v(frand(-6, 6), frand(-6, 6), frand(-6, 6));
    q->A.v = am_v(frand(-6, 6), frand(-6, 6), frand(-6, 6));
    q->B.v = am_v(frand(-6, 6), frand(-6, 6), frand(-6, 6));
    double vn = am_dot(am_sub(q->B.v, q->A.v), q->c.n);
    q->B.v = am_add(q->B.v, am_scale(q->c.n, -frand(1.0, 14.0) - vn));
}

/* unit vector perpendicular to n with random azimuth */
static am_vec3 perp_arm(am_vec3 n, double h) {
    am_vec3 ref = (fabs(n.x) < 0.9) ? am_v(1, 0, 0) : am_v(0, 1, 0);
    am_vec3 cn = am_cross(n, ref);
    am_vec3 u = am_scale(cn, 1.0 / am_len(cn));
    am_vec3 v = am_cross(n, u);
    double a = frand(0.0, 6.283185307179586);
    return am_add(am_scale(u, h * cos(a)), am_scale(v, h * sin(a)));
}

static int open_path(char* buf, size_t n, const char* dir, const char* name, FILE** fp) {
    snprintf(buf, n, "%s/%s", dir, name);
    *fp = fopen(buf, "w");
    if (!*fp) { fprintf(stderr, "cannot write %s\n", buf); return -1; }
    return 0;
}

/* precision configurations shared by the error and timing curves */
enum { K_COARSE, K_ARMERR, K_F32, K_F64, K_GUARD };
typedef struct { const char* label; int kind; double h_mm; } Cfg;

static am_solution run_cfg(Case* q, const Cfg* g) {
    double h = g->h_mm * 1e-3;
    if (g->kind == K_COARSE)
        return am_solve_apply(&q->A, &q->B, &q->c, AM_TIER_COARSE, am_v(0,0,0), am_v(0,0,0));
    if (g->kind == K_F32)
        return am_solve_apply(&q->A, &q->B, &q->c, AM_TIER_F32, am_v(0,0,0), am_v(0,0,0));
    /* K_ARMERR / K_F64 / K_GUARD all run the f64 kernel; ARMERR injects a
     * fixed arm offset (operation count identical to exact f64). */
    return am_solve_apply(&q->A, &q->B, &q->c, AM_TIER_F64,
                          g->kind == K_ARMERR ? am_v(0, h, 0) : am_v(0,0,0),
                          am_v(0,0,0));
}

int main(int argc, char** argv) {
    const char* dir = (argc > 1) ? argv[1] : ".";
    char path[4096];
    FILE* fp;
    int k;

    /* canonical calibration numbers (real, solved) */
    am_body cA0, cB0; am_contact cc0;
    am_canonical(&cA0, &cB0, &cc0);
    am_snapshot cs = am_snapshot_make(&cA0, &cB0);
    am_solution cj = am_solve_apply(&cA0, &cB0, &cc0, AM_TIER_F64, am_v(0, 0, 0), am_v(0, 0, 0));
    double j0 = cj.j, gross0 = cs.gross;
    double per_mm = j0 * 0.001 / gross0;
    double dr_max_mm = AM_TOL_ONE_PER_MILLE / per_mm;
    printf("canonical: j=%.6f N s gross=%.6f  dL_per_mm_abs=%.6e kgm2/s\n",
           j0, gross0, j0 * 0.001);
    printf("per_mm_relative=%.6e (%.4f permille/mm)  dr_max=%.4f mm\n",
           per_mm, per_mm * 1000.0, dr_max_mm);

    /* ---------------- sensitivity sweep ---------------- */
    const double hs[] = { 0.0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0 };
    const int NH = (int)(sizeof(hs) / sizeof(hs[0]));
    const int MRAND = 30000;

    if (open_path(path, sizeof path, dir, "am_sensitivity.csv", &fp)) return 2;
    fprintf(fp, "h_mm,canonical_rel,analytic_bound_rel,pair_max_rel,pair_rms_rel\n");
    int hi;
    for (hi = 0; hi < NH; ++hi) {
        double h = hs[hi], hm = h * 1e-3;
        am_body A, B; am_contact c;
        am_canonical(&A, &B, &c);
        am_snapshot s0 = am_snapshot_make(&A, &B);
        am_solution sol = am_solve_apply(&A, &B, &c, AM_TIER_F64, am_v(0, hm, 0), am_v(0, 0, 0));
        am_report rr = am_check(&s0, &A, &B, sol.eccA, sol.eccB, AM_TOL_ONE_PER_MILLE);
        double canon_rel = rr.rel_gross;
        double analytic = j0 * hm / gross0;

        double mx = 0.0, sumsq = 0.0; int cnt = 0;
        for (k = 0; k < MRAND; ++k) {
            Case q; random_case(&q);
            am_snapshot sp = am_snapshot_make(&q.A, &q.B);
            if (sp.gross < 1e-3) continue;
            am_vec3 eA = (h > 0.0) ? perp_arm(q.c.n, hm) : am_v(0, 0, 0);
            am_vec3 eB = (h > 0.0) ? perp_arm(q.c.n, hm) : am_v(0, 0, 0);
            am_solution s = am_solve_apply(&q.A, &q.B, &q.c, AM_TIER_F64, eA, eB);
            if (s.status != 1) continue;
            am_report z = am_check(&sp, &q.A, &q.B, s.eccA, s.eccB, AM_TOL_ONE_PER_MILLE);
            if (z.rel_gross > mx) mx = z.rel_gross;
            sumsq += z.rel_gross * z.rel_gross; cnt++;
        }
        double rms = cnt ? sqrt(sumsq / cnt) : 0.0;
        fprintf(fp, "%.4g,%.8e,%.8e,%.8e,%.8e\n", h, canon_rel, analytic, mx, rms);
        printf("sweep h=%5.2fmm canon=%.3e bound=%.3e pairmax=%.3e rms=%.3e n=%d\n",
               h, canon_rel, analytic, mx, rms, cnt);
    }
    fclose(fp);

    /* ---------------- shared bank for tier error + timing ---------------- */
    const int T = 4096;
    Case* bank = (Case*)malloc(sizeof(Case) * T);
    if (!bank) return 2;
    for (k = 0; k < T; ++k) random_case(&bank[k]);

    const Cfg cfgs[] = {
        { "coarse_no_cross", K_COARSE, 0.0 },
        { "arm_err_3.0mm",   K_ARMERR, 3.0 },
        { "arm_err_1.0mm",   K_ARMERR, 1.0 },
        { "arm_err_0.5mm",   K_ARMERR, 0.5 },
        { "arm_err_0.1mm",   K_ARMERR, 0.1 },
        { "float32_full",    K_F32,    0.0 },
        { "float64_full",    K_F64,    0.0 },
        { "f64_plus_guard",  K_GUARD,  0.0 }
    };
    const int NC = (int)(sizeof(cfgs) / sizeof(cfgs[0]));

    if (open_path(path, sizeof path, dir, "am_tradeoff.csv", &fp)) { free(bank); return 2; }
    fprintf(fp, "x,tier,max_rel,median_rel,passes_1e3,ns_per_collision\n");

    const int SWEEPS = 256;          /* 4096*256 ~= 1.05M collisions / round */
    const int ROUNDS = 9;

    int ci;
    for (ci = 0; ci < NC; ++ci) {
        const Cfg* g = &cfgs[ci];

        /* error over the SAME bank for every config */
        double* vals = (double*)malloc(sizeof(double) * T);
        double mx = 0.0; int nval = 0, passes = 1;
        for (k = 0; k < T; ++k) {
            Case q = bank[k];
            am_snapshot sp = am_snapshot_make(&q.A, &q.B);
            am_solution s = run_cfg(&q, g);
            if (s.status != 1 || sp.gross < 1e-3) continue;
            am_report z = am_check(&sp, &q.A, &q.B, s.eccA, s.eccB, AM_TOL_ONE_PER_MILLE);
            vals[nval++] = z.rel_gross;
            if (z.rel_gross > mx) mx = z.rel_gross;
            if (!z.pass) passes = 0;
        }
        double med = 0.0;
        if (nval > 0) { qsort(vals, nval, sizeof(double), cmp_double); med = vals[nval / 2]; }
        free(vals);

        /* timing: median of ROUNDS, fresh copy of each bank case */
        double times[ROUNDS];
        int rr;
        for (rr = 0; rr < ROUNDS; ++rr) {
            double acc = 0.0;
            double t0 = now_s();
            int sw, i;
            for (sw = 0; sw < SWEEPS; ++sw) {
                for (i = 0; i < T; ++i) {
                    Case q = bank[i];
                    double extra = 0.0;
                    am_solution s;
                    if (g->kind == K_GUARD) {
                        am_snapshot sp = am_snapshot_make(&q.A, &q.B);
                        Cfg f64 = { "f64", K_F64, 0.0 };
                        s = run_cfg(&q, &f64);
                        am_report z = am_check(&sp, &q.A, &q.B, s.eccA, s.eccB, AM_TOL_ONE_PER_MILLE);
                        extra = z.rel_gross;
                    } else {
                        s = run_cfg(&q, g);
                    }
                    acc += s.j + q.A.v.x + q.B.w.z + extra;
                }
            }
            times[rr] = (now_s() - t0) / ((double)T * SWEEPS) * 1e9;
            g_sink = acc;
        }
        /* best (least-interrupted) round is the honest throughput cost */
        double ns = times[0];
        for (rr = 1; rr < ROUNDS; ++rr) if (times[rr] < ns) ns = times[rr];

        fprintf(fp, "%d,%s,%.8e,%.8e,%d,%.3f\n", ci, g->label, mx, med, passes, ns);
        printf("tier[%d] %-18s max_rel=%.3e med=%.3e pass=%d ns/coll=%.2f\n",
               ci, g->label, mx, med, passes, ns);
    }
    fclose(fp);

    /* ---------------- reverse-solved arm-error envelope ----------------
     * Per collision: rel ~= (J/L_ref) dR, so dR_max = 1e-3 / (J/L_ref).
     * In mm this simplifies to dR_max[mm] = 1 / leverage[m^-1]. */
    if (open_path(path, sizeof path, dir, "am_design_budget.csv", &fp)) { free(bank); return 2; }
    fprintf(fp, "scope,leverage_J_over_Lref_per_m,dr_allowance_mm\n");
    {
        const int ND = 200000;
        double* lev = (double*)malloc(sizeof(double) * ND);
        int n = 0; double lmax = 0.0;
        for (k = 0; k < ND; ++k) {
            Case q; random_case(&q);
            am_snapshot sp = am_snapshot_make(&q.A, &q.B);
            if (sp.gross < 1e-3) continue;
            am_solution s = am_solve_apply(&q.A, &q.B, &q.c, AM_TIER_F64, am_v(0,0,0), am_v(0,0,0));
            if (s.status != 1 || !(s.j > 0.0)) continue;
            double z = s.j / sp.gross;
            lev[n++] = z;
            if (z > lmax) lmax = z;
        }
        qsort(lev, n, sizeof(double), cmp_double);
        double lp99 = lev[(int)(0.99 * (n - 1))];
        double lmed = lev[n / 2];
        fprintf(fp, "canonical_worked,%.6f,%.4f\n", j0 / gross0, gross0 / j0 * 1000.0 * AM_TOL_ONE_PER_MILLE);
        fprintf(fp, "random_median,%.6f,%.4f\n", lmed, 1.0 / lmed);
        fprintf(fp, "random_p99,%.6f,%.4f\n", lp99, 1.0 / lp99);
        fprintf(fp, "random_worst_sampled,%.6f,%.4f\n", lmax, 1.0 / lmax);
        fclose(fp);
        printf("envelope leverage  med=%.3f p99=%.3f worst=%.3f [1/m] -> "
               "dRmax med=%.3f p99=%.3f worst=%.4f mm\n",
               lmed, lp99, lmax, 1.0 / lmed, 1.0 / lp99, 1.0 / lmax);
        free(lev);
    }

    /* ---------------- micro timing ---------------- */
    if (open_path(path, sizeof path, dir, "am_micro_timing.csv", &fp)) { free(bank); return 2; }
    fprintf(fp, "kernel,ns_per_op\n");
    const int NX = 20000000;
    am_vec3 va[256], vb[256];
    for (k = 0; k < 256; ++k) {
        va[k] = am_v(frand(-1,1), frand(-1,1), frand(-1,1));
        vb[k] = am_v(frand(-1,1), frand(-1,1), frand(-1,1));
    }
    double acc = 0.0;
    double t0 = now_s();
    for (k = 0; k < NX; ++k) {
        am_vec3 q = am_cross(va[k & 255], vb[(k * 3 + 5) & 255]);
        acc += q.x + q.y + q.z;
    }
    double ns_cross = (now_s() - t0) / NX * 1e9;
    g_sink = acc;
    fprintf(fp, "one_am_cross,%.3f\n", ns_cross);
    printf("micro one_cross=%.3f ns\n", ns_cross);

    /* 2D A/B over runtime inputs so the compiler cannot constant-fold:
     * linear-only response vs adding the single lever-arm scalar cross. */
    const int N2 = 30000000;
    double ryv[1024];
    for (k = 0; k < 1024; ++k) ryv[k] = frand(0.05, 0.5);
    acc = 0.0;
    t0 = now_s();
    for (k = 0; k < N2; ++k) {
        double m = 1.0, e = 0.5, nx = 1.0, ny = 0.0;
        double vx = 8.0, vy = -10.0;
        double vn = vx * nx + vy * ny;
        double j = -(1 + e) * vn / (1.0 / m);
        vx += j * nx / m; vy += j * ny / m;
        acc += vx + vy + ryv[k & 1023];
    }
    double ns_lin = (now_s() - t0) / N2 * 1e9;
    g_sink = acc;

    acc = 0.0;
    t0 = now_s();
    for (k = 0; k < N2; ++k) {
        double m = 1.0, ii = 0.1, e = 0.5;
        double nx = 1.0, ny = 0.0, rx = 0.5;
        double ry = ryv[k & 1023];
        double vx = 8.0, vy = -10.0, w = 0.0;
        double arm = rx * ny - ry * nx;             /* the one lever-arm cross (z) */
        double vn = vx * nx + vy * ny + w * arm;
        double keff = 1.0 / m + arm * arm / ii;
        double j = -(1 + e) * vn / keff;
        double jx = j * nx, jy = j * ny;
        vx += jx / m; vy += jy / m;
        w += (rx * jy - ry * jx) / ii;              /* r x J angular impulse */
        acc += vx + vy + w;
    }
    double ns_arm = (now_s() - t0) / N2 * 1e9;
    g_sink = acc;
    fprintf(fp, "2d_linear_only,%.3f\n2d_with_arm_cross,%.3f\nmarginal_extra_cross,%.3f\n",
            ns_lin, ns_arm, ns_arm - ns_lin);
    printf("micro 2d_linear=%.3f ns 2d_arm=%.3f ns marginal=%.3f ns/collision\n",
           ns_lin, ns_arm, ns_arm - ns_lin);
    fclose(fp);

    free(bank);
    printf("bench_ok=1\n");
    return 0;
}
