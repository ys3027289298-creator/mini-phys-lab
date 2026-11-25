/*
 * penetration: numerical acceptance spec for "the ball must not pass through
 * the floor", evaluated at three camera distances.
 *
 * The requirement alone carries no tolerance, so the module first pins down
 * what "penetration" means (src/pen_depth.h), then derives a per-tier residual
 * quota from a screen-space visibility limit (src/pen_model.h), and finally
 * checks that the quota produces zero visible penetration frames in a
 * 60-second run at each tier.
 *
 * Constants that ARE the spec live in one block below.  Everything else is
 * derived; change a constant and the printed table and the checks move with
 * it.
 */
#include <stdio.h>
#include "pen_model.h"

/* ------------------------------------------------------- spec constants -- */
#define CAMERA_H        1.50    /* camera height above the floor, m          */
#define SENSOR_PX       1080.0  /* vertical resolution                       */
#define VFOV_DEG        60.0    /* vertical field of view                    */
#define RADIUS          0.50    /* ball radius, m                            */
#define PLANE_Y         0.0     /* floor plane, m                            */
#define GRAVITY         9.8     /* m/s^2                                     */
#define DROP            6.00    /* initial drop, m                           */
#define RESTITUTION     0.20    /* coeff. of restitution on the floor        */
#define HZ              60      /* frames per second                         */
#define K_CAP_PX        0.50    /* acceptable visual error: 0.5 px on screen */
#define K_VIS_PX        1.00    /* >= 1.0 px on screen counts as visible     */
#define R_FRAC_MAX      0.25    /* sanity guard: never sink >25% of radius   */
#define QUOTA_STEP_MM   0.0     /* mm grid for the quota; 0 = continuous     */

#define N_TIERS 3
static const double TIER_D[N_TIERS] = { 3.0, 10.0, 30.0 };
static const char*  TIER_NAME[N_TIERS] = { "3m", "10m", "30m" };

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    World w = { RADIUS, PLANE_Y, GRAVITY, CAMERA_H, SENSOR_PX, VFOV_DEG, HZ };
    double f = pen_focal_px(w.sensor_px, w.vfov_deg);
    Tier t[N_TIERS];
    for (int i = 0; i < N_TIERS; i++) {
        pen_tier(&t[i], TIER_D[i], &w, f, K_VIS_PX, K_CAP_PX);
        double allow = t[i].allow;
        if (QUOTA_STEP_MM > 0.0) {
            double step = QUOTA_STEP_MM / 1000.0;
            allow = floor(allow / step) * step;     /* round DOWN to the grid */
        }
        if (allow > R_FRAC_MAX * w.radius) allow = R_FRAC_MAX * w.radius;
        t[i].allow = allow;
    }

    /* -- 1. depth definition ------------------------------------------- */
    check("depth_positive_when_sunk", pen_depth(0.40, RADIUS, PLANE_Y) > 0.0);
    check("depth_zero_at_rest",
          fabs(pen_depth(PLANE_Y + RADIUS, RADIUS, PLANE_Y)) < 1e-12);
    check("depth_negative_when_clear",
          pen_depth(PLANE_Y + RADIUS + 0.01, RADIUS, PLANE_Y) < 0.0);

    /* -- 2. the quota is a screen fraction, not a fixed length --------- */
    check("quota_strictest_near",
          t[0].allow < t[1].allow && t[1].allow < t[2].allow);
    check("quota_below_visible_all",
          t[0].allow < t[0].A_vis && t[1].allow < t[1].A_vis &&
          t[2].allow < t[2].A_vis);
    check("quota_under_radius_ceiling",
          t[2].allow <= R_FRAC_MAX * w.radius + 1e-12);
    check("quota_scales_with_slant_distance",
          fabs(t[2].allow / t[2].R - t[0].allow / t[0].R) < 1e-12);

    /* -- 3. per-tier runs ---------------------------------------------- */
    printf("f_px=%.2f\n", f);
    long on_corr_total = 0, off_corr_total = 0, visible_total = 0;
    for (int i = 0; i < N_TIERS; i++) {
        const Tier* ti = &t[i];

        /* baseline: correction every frame (what ../integrate does today) */
        Run on = pen_run(ti, 0.0, 1.00, RESTITUTION, DROP, &w);
        /* residual policy: same integrate, deadband quota, snap correction */
        Run off = pen_run(ti, ti->allow, 1.00, RESTITUTION, DROP, &w);

        printf("tier=%s R=%.4f f=%.2f A_vis=%.3fmm quota=%.3fmm(%.3f%%) "
               "on={vis=%ld corr=%ld} off={vis=%ld corr=%ld peak=%.3fmm %.3fpx}\n",
               TIER_NAME[i], ti->R, f, ti->A_vis * 1e3, ti->allow * 1e3,
               100.0 * ti->allow / w.radius,
               on.visible_episodes, on.correction_frames,
               off.visible_episodes, off.correction_frames,
               off.max_render * 1e3, off.max_render_px);

        char nm[96];
        snprintf(nm, sizeof nm, "pen_%s_quota_zero_visible", TIER_NAME[i]);
        check(nm, off.visible_episodes == 0 && off.visible_frames == 0);

        snprintf(nm, sizeof nm, "pen_%s_render_under_quota", TIER_NAME[i]);
        check(nm, off.max_render <= ti->allow + 1e-12);

        snprintf(nm, sizeof nm, "pen_%s_render_under_one_pixel", TIER_NAME[i]);
        check(nm, off.max_render_px <= K_VIS_PX);

        /* curve trade-off: stretching the quota must not create corrections */
        snprintf(nm, sizeof nm, "pen_%s_not_more_corrections", TIER_NAME[i]);
        check(nm, off.correction_frames <= on.correction_frames);

        /* negative control: the measurement must SEE a sink when the quota is
         * deliberately set past the visible limit.  Without this the zero
         * above could be a test that cannot fail. */
        Run bad = pen_run(ti, 3.0 * ti->A_vis, 1.00, RESTITUTION, DROP, &w);
        snprintf(nm, sizeof nm, "pen_%s_detects_overcap_sink", TIER_NAME[i]);
        check(nm, bad.visible_episodes > 0);

        /* preset values from the derivation must reproduce exactly */
        visible_total += off.visible_frames;
        on_corr_total += on.correction_frames;
        off_corr_total += off.correction_frames;
    }

    /* -- 4. residual vs correction coefficient ------------------------- */
    /*
     * Per frame the integrator adds one free-fall increment
     *     d = g*dt^2
     * to the depth (velocity is cleared on contact, so gravity re-adds d
     * every frame).  A correction that removes a fraction beta of the
     * measured depth h leaves
     *     h_next = h - beta*h + d  =  (1-beta)*h + d.
     * The fixed point is h* = d/beta, so the residual that survives a
     * correction frame is
     *
     *     residual = (1-beta)*h* = ((1-beta)/beta) * d.
     *
     * beta = 1 -> residual 0 (the corrected frame is clean).  beta < 1 is
     * strictly worse and diverges as beta -> 0: beta = 1/2 leaves d,
     * beta = 1/4 leaves 3d.  So the spec must use beta = 1, and the only
     * knob left is the deadband quota.  Verify the closed form against the
     * simulation's settled frame.
     */
    double dt = 1.0 / HZ;
    double d = GRAVITY * dt * dt;
    double betas[3] = { 1.00, 0.50, 0.25 };
    for (int k = 0; k < 3; k++) {
        Run rr = pen_run(&t[0], 0.0, betas[k], 0.0, d, &w);
        double want = ((1.0 - betas[k]) / betas[k]) * d;
        char nm[96];
        snprintf(nm, sizeof nm, "residual_matches_coefficient_beta_%d",
                 (int)(betas[k] * 100));
        check(nm, fabs(rr.final_render - want) < 1e-9);
    }

    /* -- 5. residual <-> correction frequency -------------------------- */
    /*
     * The deadband policy corrects only once the depth has grown past the
     * quota `allow`.  Starting from 0 after a correction, depth grows by
     * d = g*dt^2 per frame (velocity is cleared on contact), so a correction
     * fires every
     *
     *     m = ceil(allow / d)                       frames
     *
     * and the fraction of frames touched is 1/m.  Two curves therefore move
     * in opposite directions as `allow` grows:
     *
     *     visible(allow)    0 while allow < A_vis, then jumps at A_vis
     *     cost(allow)       1/m = d/allow, falling like 1/allow
     *
     * The chosen quota is the largest allow that keeps visible() at zero,
     * i.e. just under A_vis, which also minimises cost.  Check the two against
     * the simulation and print them.
     */
    double dt5 = 1.0 / HZ;
    double d5 = GRAVITY * dt5 * dt5;
    printf("frame_cost: depth_per_frame=%.4fmm\n", d5 * 1e3);
    for (int i = 0; i < N_TIERS; i++) {
        long m = (long)ceil(t[i].allow / d5);
        if (m < 1) m = 1;
        /* Sustained-contact regime: e = 0 so the ball stays on the floor
         * between corrections.  The cost is a *rate* (fraction of frames
         * corrected); the absolute count also absorbs the initial drop. */
        Run rr = pen_run(&t[i], t[i].allow, 1.00, 0.0, DROP, &w);
        double frac = (double)rr.correction_frames / (double)rr.frames;
        double want = 1.0 / (double)m;
        printf("cost_curve tier=%s m=%ld want_frac=%.4f measured_frac=%.4f "
               "corrections=%ld/%ld\n",
               TIER_NAME[i], m, want, frac, rr.correction_frames, rr.frames);
        char nm[96];
        snprintf(nm, sizeof nm, "pen_%s_cost_matches_formula", TIER_NAME[i]);
        check(nm, fabs(frac - want) <= 0.02);
    }

    /* -- 6. frame-rate robustness -------------------------------------- */
    /*
     * With beta = 1 the rendered depth never exceeds the quota regardless of
     * d, so visibility is a function of quota and distance only.  Only the
     * correction cost moves with frame rate.  Verify at 30 and 120 Hz.
     */
    int hz_list[3] = { 30, 60, 120 };
    for (int j = 0; j < 3; j++) {
        World wj = w;
        wj.hz = hz_list[j];
        for (int i = 0; i < N_TIERS; i++) {
            Run rj = pen_run(&t[i], t[i].allow, 1.00, RESTITUTION, DROP, &wj);
            char nm[96];
            snprintf(nm, sizeof nm, "pen_%s_%dhz_zero_visible", TIER_NAME[i], hz_list[j]);
            check(nm, rj.visible_frames == 0);
        }
    }

    /* -- 7. scope of the guarantee ------------------------------------- */
    /*
     * The quota bounds what is RENDERED: correction runs before the frame is
     * presented, so a viewer never sees more than `allow` whether the ball
     * just landed or has been resting.  It does NOT bound the pre-correction
     * depth: at 10.8 m/s impact the ball travels ~181 mm per 60 Hz frame, so
     * the raw depth on the landing frame is anywhere in [0, 181] mm depending
     * on where the frame boundary falls.  That is still < one diameter
     * (1000 mm), so the ball never ends up on the far side of the floor, but
     * the number is frame-phase dependent and must not be quoted as a spec.
     */
    {
        Run land = pen_run(&t[0], 0.0, 1.00, RESTITUTION, DROP, &w);
        double v_impact = sqrt(2.0 * GRAVITY * DROP);
        double travel = v_impact / HZ;
        printf("impact: v=%.2fm/s per_frame_travel=%.1fmm diameter=%.1fmm "
               "raw_physical_max=%.2fmm rendered_max=%.3fmm\n",
               v_impact, travel * 1e3, 2.0 * w.radius * 1e3,
               land.max_physical * 1e3, land.max_render * 1e3);
        check("physical_stays_under_diameter", land.max_physical < 2.0 * w.radius);
        check("rendered_is_bounded_at_impact", land.max_render < 1e-12);
    }

    /* -- 8. reported quota values -------------------------------------- */
    printf("spec_quota_mm: 3m=%.3f 10m=%.3f 30m=%.3f\n",
           t[0].allow * 1e3, t[1].allow * 1e3, t[2].allow * 1e3);
    printf("tot_visible=%ld tot_corr_on=%ld tot_corr_off=%ld\n",
           visible_total, on_corr_total, off_corr_total);
    check("total_visible_zero", visible_total == 0);
    check("residual_never_costs_more", off_corr_total <= on_corr_total);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
