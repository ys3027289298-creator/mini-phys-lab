
/* integrate: semi-implicit euler with position correction against penetration.
 *
 * Substep window for a thin floor:
 *
 *     anti-tunnelling:   s_max = V*T/N < h
 *     anti-saturation:   g*(T/N)^2 > eps       (zero-speed/rest substep)
 *
 * The first inequality gives a lower bound on N; the second gives an upper
 * bound.  substep_window() contains the constant-speed form, while
 * gravity_live_limit() contains the acceleration-aware form.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

typedef struct { float x, y, px, py, vx, vy; } Body;

static int step(Body* b, float dt, float gx, float gy, float floor_y, float radius) {
    if (dt < 0.0f || radius <= 0.0f) return -1;
    b->px = b->x; b->py = b->y;
    b->vx += gx * dt;
    b->vy += gy * dt;
    b->x += b->vx * dt;
    b->y += b->vy * dt;
    /* floor at y = floor_y, body bottom is y - radius */
    float pen = (floor_y + radius) - b->y;
    if (pen > 0.0f) {
        b->y += pen;              /* positional correction */
        if (b->vy < 0.0f) b->vy = 0.0f;
    }
    return pen > 0.0f ? 1 : 0;
}

/* ---------------- thin-floor substep study ---------------- */

#define STUDY_G        9.81f
#define STUDY_Y0       5.0f
#define STUDY_H        0.10f
#define STUDY_R        0.02f
#define STUDY_EPS      1.0e-4f
#define STUDY_RUN_S    1.20f
#define STUDY_LIMIT    0.05f
#define STUDY_SWEEP_N  40

typedef struct {
    float top_y;      /* y of the upper face */
    float thickness;  /* slab thickness h */
    float radius;     /* ball radius r */
    float eps;        /* displacements with |dy| < eps are snapped to zero */
    float gravity_y;  /* signed gravity acceleration */
} ThinFloor;

typedef struct {
    int substeps;
    int zero_steps;       /* actual position increment was zero */
    int truncated_steps;  /* raw nonzero increment was snapped to zero */
    int tunnel_count;     /* one step jumped completely over the detector band */
    int contact_steps;
    float max_step;
} FrameStats;

/* Constant-speed integer window.
 *
 * V*T/N < h   -> N >= floor(V*T/h) + 1
 * V*T/N > eps -> N <= ceil(V*T/eps) - 1
 *
 * Returns 0 when an integer window exists, 1 when it is empty, <0 on bad args.
 */
static int substep_window(double frame_dt, double max_speed,
                          double h, double eps,
                          int* n_safe, int* n_speed_live) {
    if (!(frame_dt > 0.0 && max_speed >= 0.0 && h > 0.0 && eps > 0.0)) return -1;

    int lo = (int)floor(max_speed * frame_dt / h) + 1;
    int hi = (int)ceil(max_speed * frame_dt / eps) - 1;
    if (lo < 1) lo = 1;
    if (n_safe) *n_safe = lo;
    if (n_speed_live) *n_speed_live = hi;
    return lo <= hi ? 0 : 1;
}

/* Acceleration-aware upper bound.  From zero vertical speed, semi-implicit
 * Euler's first displacement is g*tau^2.  Requiring it to exceed eps gives
 *
 *     N < T*sqrt(g/eps),   N_max = ceil(T*sqrt(g/eps)) - 1.
 */
static int gravity_live_limit(double frame_dt, double g, double eps, int* n_max) {
    if (!(frame_dt > 0.0 && g > 0.0 && eps > 0.0)) return -1;
    double root = frame_dt * sqrt(g / eps);
    if (n_max) *n_max = (int)ceil(root) - 1;
    return 0;
}

/* Exact endpoint-skip geometry for a sphere against slab [top-h, top].
 * The centre detects overlap on [top-h-r, top+r], width h+2r.  A downward
 * step of length s can have both endpoints outside iff s > h+2r.
 */
static int endpoint_skip_possible(double step_len, double h, double r) {
    return step_len > h + 2.0 * r;
}

/* Swept segment vs the expanded centre band [lo, hi].  This is the ray/CCD
 * alternative: one large frame step is still caught by intersecting the whole
 * segment rather than testing only its endpoints.
 */
static int swept_interval_hit(double y0, double y1, double lo, double hi,
                              double* enter_t) {
    double dy = y1 - y0;
    if (fabs(dy) < 1.0e-15) {
        int hit = y0 >= lo && y0 <= hi;
        if (hit && enter_t) *enter_t = 0.0;
        return hit;
    }

    double ta = (lo - y0) / dy;
    double tb = (hi - y0) / dy;
    double t_lo = fmin(ta, tb);
    double t_hi = fmax(ta, tb);
    double enter = fmax(t_lo, 0.0);
    double leave = fmin(t_hi, 1.0);
    if (enter <= leave + 1.0e-12) {
        if (enter_t) *enter_t = enter;
        return 1;
    }
    return 0;
}

static int thin_floor_step(Body* b, float tau, const ThinFloor* f, FrameStats* st) {
    if (tau <= 0.0f || f->thickness <= 0.0f || f->radius <= 0.0f || f->eps <= 0.0f)
        return -1;

    b->px = b->x;
    b->py = b->y;
    b->vy += f->gravity_y * tau;

    float raw_dy = b->vy * tau;
    float dy = fabsf(raw_dy) < f->eps ? 0.0f : raw_dy;
    if (dy == 0.0f) {
        st->zero_steps++;
        if (raw_dy != 0.0f) st->truncated_steps++;
    }

    b->y += dy;
    st->substeps++;
    float moved = fabsf(dy);
    if (moved > st->max_step) st->max_step = moved;

    double upper = (double)f->top_y + f->radius;
    double lower = (double)f->top_y - f->thickness - f->radius;

    int crossed = (double)b->py > upper && (double)b->y < lower;
    if (crossed) st->tunnel_count++;

    int hit = (double)b->y >= lower && (double)b->y <= upper;
    if (hit) {
        b->y = f->top_y + f->radius;
        if (b->vy < 0.0f) b->vy = 0.0f;
        st->contact_steps++;
    }
    return crossed ? 1 : (hit ? 2 : 0);
}

static int step_substepped(Body* b, float frame_dt, int n,
                           const ThinFloor* f, FrameStats* st) {
    if (frame_dt <= 0.0f || n <= 0) return -1;
    float tau = frame_dt / (float)n;
    for (int i = 0; i < n; i++) {
        int rc = thin_floor_step(b, tau, f, st);
        if (rc < 0) return rc;
    }
    return 0;
}

static int run_study_tier(float frame_dt, int n, FrameStats* st) {
    memset(st, 0, sizeof(*st));
    Body b;
    memset(&b, 0, sizeof(b));
    b.y = STUDY_Y0;

    ThinFloor floor = {
        0.0f, STUDY_H, STUDY_R, STUDY_EPS, -STUDY_G
    };
    int frames = (int)lroundf(STUDY_RUN_S / frame_dt);
    for (int frame = 0; frame < frames; frame++) {
        int rc = step_substepped(&b, frame_dt, n, &floor, st);
        if (rc != 0) return rc;
    }
    return 0;
}

/* Worst-phase discrete miss count, mirroring tools/substep_study.py. */
static int phase_miss_count(double frame_dt, int n, double speed,
                            double h, double r, int phases) {
    double s = speed * frame_dt / (double)n;
    double upper = r;
    double lower = -h - r;
    int misses = 0;

    for (int j = 0; j < phases; j++) {
        double phase = ((double)j + 0.5) / (double)phases * s;
        double y_before = upper + phase;
        double y_after = y_before - s;
        int hit_before = y_before >= lower && y_before <= upper;
        int hit_after = y_after >= lower && y_after <= upper;
        if (!hit_before && !hit_after && y_after < lower) misses++;
    }
    return misses;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Body b; memset(&b, 0, sizeof(b));
    b.y = 10.0f;
    int hit = step(&b, 0.1f, 0, -9.8f, 0.0f, 0.5f);
    check("falling_no_hit", hit == 0 && b.y < 10.0f && b.y > 0.0f);

    /* drive it into the floor */
    for (int i = 0; i < 100; i++) step(&b, 0.1f, 0, -9.8f, 0.0f, 0.5f);
    check("lands_on_floor", fabsf(b.y - 0.5f) < 1e-3f);
    check("vertical_stopped", fabsf(b.vy) < 1e-6f);

    /* never sinks below floor across a long run */
    Body c; memset(&c, 0, sizeof(c));
    c.y = 3.0f; c.vx = 5.0f;
    float min_y = 1e9f;
    for (int i = 0; i < 500; i++) {
        step(&c, 0.05f, 0, -9.8f, 0.0f, 0.5f);
        if (c.y < min_y) min_y = c.y;
    }
    check("never_sinks", min_y >= 0.5f - 1e-4f);

    Body d; memset(&d, 0, sizeof(d));
    check("bad_dt", step(&d, -1.0f, 0, 0, 0, 0.5f) == -1);
    check("bad_radius", step(&d, 0.1f, 0, 0, 0, 0.0f) == -1);

    /* ---------------- substep derivation and three tiers ---------------- */

    double v_impact = sqrt(2.0 * (double)STUDY_G * (double)STUDY_Y0);
    float dts[3] = { 1.0f / 60.0f, 1.0f / 30.0f, 0.10f };
    int expected_n[3] = { 2, 4, 10 };
    int expected_gmax[3] = { 5, 10, 31 };
    int probe_n[3] = { 8, 16, 32 };

    printf("\nsubstep study: V=%.4f h=%.2f r=%.2f eps=%.0e\n",
           v_impact, STUDY_H, STUDY_R, STUDY_EPS);
    printf("strict tau window: %.3f ms < tau < %.3f ms\n",
           1000.0 * sqrt((double)STUDY_EPS / (double)STUDY_G),
           1000.0 * (double)STUDY_H / v_impact);

    int all_tiers_ok = 1;
    for (int tier = 0; tier < 3; tier++) {
        float dt = dts[tier];
        int n_min = 0, n_speed_max = 0, n_gmax = 0;
        int window_rc = substep_window(dt, v_impact, STUDY_H, STUDY_EPS,
                                       &n_min, &n_speed_max);
        gravity_live_limit(dt, STUDY_G, STUDY_EPS, &n_gmax);

        int n5_scan = 0;
        for (int n = 1; n <= STUDY_SWEEP_N; n++) {
            FrameStats scanned;
            run_study_tier(dt, n, &scanned);
            float ratio = (float)scanned.zero_steps / (float)scanned.substeps;
            if (ratio <= STUDY_LIMIT) n5_scan = n;
        }

        FrameStats chosen;
        run_study_tier(dt, expected_n[tier], &chosen);
        int phase_miss = phase_miss_count(dt, expected_n[tier], v_impact,
                                          STUDY_H, STUDY_R, 401);
        float ratio = (float)chosen.zero_steps / (float)chosen.substeps;
        int tier_ok = window_rc == 0 &&
                      n_min == expected_n[tier] &&
                      n_gmax == expected_gmax[tier] &&
                      n5_scan == expected_gmax[tier] &&
                      phase_miss == 0 &&
                      chosen.tunnel_count == 0 &&
                      ratio <= STUDY_LIMIT;
        all_tiers_ok &= tier_ok;

        printf("T_ms=%6.2f N=%2d tau_ms=%6.3f bounds=[%d..%d] "
               "scan_N5=%2d phase_miss=%3d tunnels=%d zero=%3d/%3d "
               "ratio=%5.2f%% calls/s=%6.1f %s\n",
               1000.0f * dt, expected_n[tier], 1000.0f * dt / expected_n[tier],
               n_min, n_gmax, n5_scan, phase_miss, chosen.tunnel_count,
               chosen.zero_steps, chosen.substeps, 100.0f * ratio,
               (float)expected_n[tier] / dt, tier_ok ? "PASS" : "FAIL");

        FrameStats probe;
        run_study_tier(dt, probe_n[tier], &probe);
        float probe_ratio = (float)probe.zero_steps / (float)probe.substeps;
        int probe_bad = probe_ratio > STUDY_LIMIT;
        printf("  over-split N=%2d: tunnels=%d zero_ratio=%5.2f%% -> %s N<=%d\n",
               probe_n[tier], probe.tunnel_count, 100.0f * probe_ratio,
               probe_bad ? "reduce to" : "within limit; no correction needed,",
               expected_gmax[tier]);
        if (!probe_bad) all_tiers_ok = 0;
    }
    check("three_frame_tiers_zero_penetration_and_le_5pct", all_tiers_ok);

    /* The constant-speed upper bound is intentionally loose; the gravity/rest
     * bound is what makes the 5% resting-tail cliff at these parameters. */
    int lo = 0, hi = 0;
    substep_window(1.0 / 60.0, v_impact, STUDY_H, STUDY_EPS, &lo, &hi);
    check("speed_bound_is_loose", lo == 2 && hi == 1650);

    /* h <= eps makes s<h and s>eps mutually impossible. */
    int empty_h_equals_eps = substep_window(0.1, v_impact, 0.001, 0.001,
                                            NULL, NULL);
    check("empty_window_when_h_le_eps", empty_h_equals_eps == 1);

    /* At H=5 m, V^2=2gH, so eps >= g*h^2/V^2 = h^2/(2H) = 1 mm closes the
     * acceleration-aware window by one strict integer step for every tier. */
    int gm_100ms_1mm = 0;
    gravity_live_limit(0.1, 9.81, 0.001, &gm_100ms_1mm);
    check("one_mm_quantization_closes_window", gm_100ms_1mm == 9);

    /* Radius geometry: r>h alone is not enough.  Both endpoint misses require
     * s > h+2r, equivalently r < (s-h)/2 for a fixed s.  With r=0.20 > h=0.10
     * and a 1.2 m step the endpoints can miss; N=10 shrinks the step below the
     * expanded 0.50 m band and ordinary substepping catches it. */
    float big_r = 0.20f;
    check("radius_gt_thickness_one_jump_can_skip",
          endpoint_skip_possible(1.2, STUDY_H, big_r) && big_r > STUDY_H);
    check("normal_substeps_catch_large_radius",
          !endpoint_skip_possible(v_impact * 0.1 / 10.0, STUDY_H, big_r));

    double enter = -1.0;
    int ray_hit = swept_interval_hit(0.21, -0.99,
                                     -STUDY_H - big_r, big_r, &enter);
    check("ray_detects_large_radius_single_jump", ray_hit && enter >= 0.0 && enter <= 1.0);

    /* If truncation is so coarse that eps >= h+2r, even the expanded band can
     * be skipped by the smallest representable move; ray on raw motion or
     * a thicker query slab is then required. */
    int exact_band_window = substep_window(0.1, v_impact,
                                           STUDY_H + 2.0 * big_r, 0.60,
                                           NULL, NULL);
    check("coarse_quantization_can_empty_expanded_band", exact_band_window == 1);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
