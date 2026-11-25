/*
 * deform: self-tests, stiffness calibration from video squash ratios and
 * three-drop-height validation.  `make test` runs everything.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "deform.h"
#include "video_measurements.h"

/* ---- shared scenario ------------------------------------------------- */

#define BALL_MASS    1.0f     /* kg  */
#define BALL_RADIUS  0.1f     /* m   */
#define GRAVITY      9.8f     /* m/s^2 */
#define RESTITUTION  0.85f    /* target bounce: h_bounce = e^2 h_drop */
#define DT           2e-5f    /* s, << contact oscillation period */
#define MAX_STEPS    2000000

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static int feq(float a, float b, float tol) {
    return fabsf(a - b) <= tol;
}

static void make_body(DeformBody* b, float k, float e) {
    b->mass = BALL_MASS;
    b->radius = BALL_RADIUS;
    b->stiffness = k;
    b->damping = deform_damping_for(BALL_MASS, k, e);
    b->gravity = GRAVITY;
}

/* ---- unit tests ------------------------------------------------------ */

static void unit_tests(void) {
    DeformBody b;
    make_body(&b, 20000.0f, RESTITUTION);

    /* validation */
    DeformBody bad = b; bad.mass = 0.0f;
    check("bad_body_rejected", deform_x_max(4.0f, &bad) < 0.0f);
    check("bad_speed_rejected", deform_x_max(0.0f, &b) < 0.0f);
    check("bad_calib_inputs",
          deform_stiffness_for(4.0f, -1.0f, 1.0f, 0.1f, 0.05f) < 0.0f);

    /* undamped energy balance: x_max = v sqrt(m/k) */
    DeformBody lossless = b; lossless.damping = 0.0f;
    float x_energy = 4.0f * sqrtf(lossless.mass / lossless.stiffness);
    check("undamped_energy_balance",
          feq(deform_x_max(4.0f, &lossless), x_energy, 1e-5f));
    check("peak_factor_undamped_one", feq(deform_peak_factor(0.0f), 1.0f, 1e-6f));

    /* restitution <-> damping ratio round trip */
    float zeta = deform_zeta_from_restitution(RESTITUTION);
    check("restitution_roundtrip",
          feq(deform_restitution(zeta), RESTITUTION, 1e-4f));
    check("damping_matches_zeta", feq(deform_zeta(&b), zeta, 1e-5f));

    /* analytic squash amplitude vs time-stepped simulation */
    ContactResult r; memset(&r, 0, sizeof(r));
    float v_mid = sqrtf(2.0f * GRAVITY * 1.0f);
    int rc = deform_simulate(&b, v_mid, DT, MAX_STEPS, &r);
    float x_analytic = deform_x_max(v_mid, &b);
    check("simulate_ok", rc == 0 && r.contact_time > 0.0f);
    check("sim_matches_analytic",
          feq(r.x_max, x_analytic, 0.03f * x_analytic));

    /* restitution observed in the simulation matches the material target;
     * the analytic e assumes zero gravity, but gravity keeps acting over
     * the ~17 ms contact, so a few-1e-3 shortfall is physical */
    check("sim_restitution",
          feq(r.exit_speed / v_mid, RESTITUTION, 5e-3f));
    check("sim_bounce_height",
          feq(r.bounce_height, RESTITUTION * RESTITUTION * 1.0f,
              1e-2f * RESTITUTION * RESTITUTION * 1.0f));

    /* squash ratio definition: x over resting diameter */
    check("squash_ratio_definition",
          feq(r.squash_ratio, r.x_max / (2.0f * BALL_RADIUS), 1e-6f));

    /* calibration inverse: k from a target ratio reproduces that ratio */
    float target = 0.113f;
    float k_cal = deform_stiffness_for(v_mid, target,
                                       BALL_MASS, BALL_RADIUS, zeta);
    DeformBody cal; make_body(&cal, k_cal, RESTITUTION);
    float x_cal = deform_x_max(v_mid, &cal);
    check("calibration_inverse",
          feq(x_cal / (2.0f * BALL_RADIUS), target, 1e-4f));

    /* volume-preserving visual scales */
    float sx, sy;
    deform_scales(0.2f, &sx, &sy);
    check("visual_scales_preserve_area", feq(sx * sx * sy, 1.0f, 1e-5f));
    deform_scales(0.0f, &sx, &sy);
    check("visual_scales_rest_round", feq(sx, 1.0f, 1e-6f) && feq(sy, 1.0f, 1e-6f));

    /* stiffer ball squashes less at the same impact */
    DeformBody stiff; make_body(&stiff, 4.0f * b.stiffness, RESTITUTION);
    check("stiffer_squashes_less",
          deform_x_max(v_mid, &stiff) < deform_x_max(v_mid, &b));

    /* over-damped contact (e ~ 0) reports no resolved bounce */
    DeformBody mushy; make_body(&mushy, 5000.0f, 0.0f);
    ContactResult rm; memset(&rm, 0, sizeof(rm));
    deform_simulate(&mushy, v_mid, DT, MAX_STEPS, &rm);
    check("overdamped_no_bounce", rm.contact_time == 0.0f && rm.bounce_height == 0.0f);
}

/* ---- calibration + validation ---------------------------------------- */

/*
 * Pick k off the middle drop using its video squash ratio, then run all
 * three drop heights.  Returns the number of tolerance failures.
 */
static int calibrate_and_validate(void) {
    float zeta = deform_zeta_from_restitution(RESTITUTION);
    const VideoMeasure* mid = 0;
    for (int i = 0; i < VIDEO_N_MEASURES; i++)
        if (VideoMeasurements[i].drop_height == 1.0f) mid = &VideoMeasurements[i];
    if (mid == 0) { printf("calib_error=middle_drop_missing\n"); return 1; }

    float v_mid = sqrtf(2.0f * GRAVITY * mid->drop_height);
    float k = deform_stiffness_for(v_mid, mid->squash_ratio,
                                   BALL_MASS, BALL_RADIUS, zeta);
    DeformBody body; make_body(&body, k, RESTITUTION);

    printf("\n--- calibration ---\n");
    printf("video_source=%s\n", VIDEO_SOURCE);
    printf("calib_drop=%.2fm v_impact=%.4fm/s video_epsilon=%.4f\n",
           mid->drop_height, v_mid, mid->squash_ratio);
    printf("zeta=%.5f restitution=%.3f\n", zeta, RESTITUTION);
    printf("stiffness_k=%.2f N/m  damping_c=%.4f Ns/m\n",
           k, body.damping);
    printf("wn=%.2f rad/s contact_period~%.4fs dt=%.1es\n",
           deform_wn(&body), (float)M_PI / deform_wn(&body), DT);

    printf("\n--- three drop heights ---\n");
    printf("%-6s %-9s %-10s %-10s %-10s %-10s %-8s %s\n",
           "h[m]", "v_imp", "eps_sim", "eps_video", "rel_dev",
           "h_bounce", "tol10%", "status");

    int n_fail = 0;
    for (int i = 0; i < VIDEO_N_MEASURES; i++) {
        const VideoMeasure* m = &VideoMeasurements[i];
        float v = sqrtf(2.0f * GRAVITY * m->drop_height);
        ContactResult r; memset(&r, 0, sizeof(r));
        deform_simulate(&body, v, DT, MAX_STEPS, &r);

        float dev = (r.squash_ratio - m->squash_ratio) / m->squash_ratio;
        int ok = fabsf(dev) <= VIDEO_TOLERANCE;
        if (!ok) n_fail++;

        printf("%-6.2f %-9.4f %-10.4f %-10.4f %-+10.4f %-10.4f %-8s %s",
               m->drop_height, v, r.squash_ratio, m->squash_ratio, dev,
               r.bounce_height, ok ? "pass" : "FAIL",
               ok ? "ok\n" : "out_of_tolerance\n");
        if (!ok) {
            /* required disclosure: actual simulated vs video-measured */
            printf("  >>> MISMATCH h=%.2fm: actual_squash=%.4f "
                   "video_squash=%.4f deviation=%.1f%% (limit %.0f%%)\n",
                   m->drop_height, r.squash_ratio, m->squash_ratio,
                   dev * 100.0f, VIDEO_TOLERANCE * 100.0f);
        }
    }
    printf("tolerance_failures=%d\n", n_fail);
    return n_fail;
}

/* CSV emitter for the two trade-off curves (swept k, middle drop). */
static int sweep_csv(const char* path) {
    FILE* f = fopen(path, "w");
    if (!f) { perror(path); return 1; }
    float v = sqrtf(2.0f * GRAVITY * 1.0f);
    fprintf(f, "k_Npm,squash_ratio,bounce_height_m\n");
    for (int i = 0; i <= 40; i++) {
        float k = powf(10.0f, 3.0f + 3.0f * i / 40.0f);   /* 1e3 .. 1e6 */
        DeformBody body; make_body(&body, k, RESTITUTION);
        ContactResult r; memset(&r, 0, sizeof(r));
        deform_simulate(&body, v, DT, MAX_STEPS, &r);
        float x = deform_x_max(v, &body);
        float eps = (x > 0.0f) ? x / (2.0f * BALL_RADIUS) : r.squash_ratio;
        fprintf(f, "%.3f,%.6f,%.6f\n", k, eps, r.bounce_height);
    }
    fclose(f);
    printf("wrote %s\n", path);
    return 0;
}

int main(int argc, char** argv) {
    if (argc == 3 && strcmp(argv[1], "--sweep") == 0)
        return sweep_csv(argv[2]);

    unit_tests();
    int n_fail = calibrate_and_validate();

    printf("\nall_ok=%d\n", failures == 0 && n_fail == 0);
    return (failures == 0 && n_fail == 0) ? 0 : 1;
}
