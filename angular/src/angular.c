
/* angular self-test + conservation guard.
 * `make test` builds this and runs it. Every collision is checked against
 * the one-per-mille inequality; failures print AM_VIOLATION with the angular
 * momentum change and the eccentricity at that collision. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "am.h"

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static void emit(const char* tag, const am_report* r) {
    if (r->pass) {
        printf("AM_PASS case=%-22s dL=%.3e rel_gross=%.3e rel_total=%.3e "
               "tol=%.0e eccA=%.4f eccB=%.4f\n",
               tag, r->dL_abs, r->rel_gross,
               isfinite(r->rel_total) ? r->rel_total : -1.0,
               r->tol, r->eccA, r->eccB);
    } else {
        /* Required alarm: report the change amount and the eccentricity. */
        printf("AM_VIOLATION case=%-18s dL_kgm2/s=%.6g rel=%.6g limit=%.0e "
               "excess_x=%.2f ecc_max_m=%.5f eccA_m=%.5f eccB_m=%.5f "
               "near_cancel=%d\n",
               tag, r->dL_abs, r->rel_gross, r->tol,
               r->rel_gross / r->tol, r->ecc_max, r->eccA, r->eccB,
               r->near_cancel);
    }
}

/* snapshot -> solve/apply -> independent check, all in one call */
static am_report run_case(const char* tag, am_tier tier,
                          am_vec3 errA, am_vec3 errB, int loud) {
    am_body A0, B0;
    am_contact c;
    am_canonical(&A0, &B0, &c);

    am_snapshot snap = am_snapshot_make(&A0, &B0);
    am_solution sol = am_solve_apply(&A0, &B0, &c, tier, errA, errB);
    am_report r = am_check(&snap, &A0, &B0, sol.eccA, sol.eccB, AM_TOL_ONE_PER_MILLE);
    if (loud) emit(tag, &r);
    return r;
}

/* deterministic xorshift for the broad battery */
static unsigned long long rng_state = 0x123456789abcdef0ULL;
static double urand(void) {
    unsigned long long x = rng_state;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    rng_state = x;
    return (double)((x >> 11) & ((1ULL << 53) - 1)) / (double)(1ULL << 53);
}
static double frand(double lo, double hi) { return lo + (hi - lo) * urand(); }

int main(void) {
    /* ---- sign convention primitives ---- */
    am_vec3 ex = am_v(1, 0, 0), ey = am_v(0, 1, 0);
    check("cross_right_hand_z", fabs(am_cross(ex, ey).z - 1.0) < 1e-15);
    check("cross_anti_commutes", fabs(am_cross(ey, ex).z + 1.0) < 1e-15);
    check("tol_is_one_per_mille", AM_TOL_ONE_PER_MILLE == 1e-3);

    /* ---- canonical hand-computed values ---- */
    am_body A, B; am_contact c;
    am_canonical(&A, &B, &c);
    am_snapshot snap0 = am_snapshot_make(&A, &B);
    check("canonical_origin_com", fabs(snap0.origin.x - 1.0 / 6.0) < 1e-12
          && fabs(snap0.origin.y) < 1e-12);
    check("canonical_Lbefore_3", fabs(snap0.L_before.z - 3.0) < 1e-12
          && snap0.total_before > 0.0);
    check("canonical_gross_3", fabs(snap0.gross - 3.0) < 1e-12);

    am_solution sol0 = am_solve_apply(&A, &B, &c, AM_TIER_F64,
                                      am_v(0, 0, 0), am_v(0, 0, 0));
    check("canonical_approaching", sol0.status == 1 && sol0.vn < 0.0);
    check("canonical_ecc", fabs(sol0.eccA - 0.30) < 1e-12
          && fabs(sol0.eccB - 0.15) < 1e-12);
    check("canonical_keff", fabs(sol0.k_eff - 2.5125) < 1e-9);
    check("canonical_j", fabs(sol0.j - 400.0 / 67.0) < 1e-9);

    /* ---- the requirement itself: exact double solve conserves < 1e-3 ---- */
    am_report good = am_check(&snap0, &A, &B, sol0.eccA, sol0.eccB, AM_TOL_ONE_PER_MILLE);
    emit("canonical_exact_f64", &good);
    check("exact_f64_conserves_1e3", good.pass && good.rel_gross < 1e-9);

    /* ---- coarse tier: full particle solve, lever arm ignored.
     * K = 1/mA+1/mB = 1.5 -> j = 15/1.5 = 10 (larger than the correct
     * 400/67); orbital change -0.45*10 = -4.5 -> rel = 1.5 (150%). ---- */
    am_report coarse = run_case("particle_no_arm", AM_TIER_COARSE,
                                am_v(0, 0, 0), am_v(0, 0, 0), 1);
    check("coarse_violates_budget", !coarse.pass);
    check("coarse_dL_matches_hand", fabs(coarse.dL_abs - 4.5) < 1e-9);
    check("coarse_reports_ecc", fabs(coarse.ecc_max - 0.30) < 1e-9);
    check("coarse_detected_as_150pct", fabs(coarse.rel_gross - 1.5) < 1e-9);

    /* ---- isolated missing r x J cross: correct j = 400/67, linear impulse
     * applied, but the spin angular impulse is omitted. Orbital change is
     * -0.45*j = -180/67 ~= -2.6866 with no compensating spin. ---- */
    {
        am_body P, Q; am_contact c2;
        am_canonical(&P, &Q, &c2);
        am_snapshot sp = am_snapshot_make(&P, &Q);
        double jcorrect = 400.0 / 67.0;
        P.v.x += -jcorrect / P.m;          /* impulse -j n on A, linear only */
        Q.v.x += +jcorrect / Q.m;          /* impulse +j n on B             */
        /* P.w / Q.w deliberately left zero: the r x J cross was skipped */
        am_report miss = am_check(&sp, &P, &Q, 0.30, 0.15, AM_TOL_ONE_PER_MILLE);
        emit("missing_rJ_crossproduct", &miss);
        check("missing_spin_cross_caught", !miss.pass);
        check("missing_spin_cross_dL", fabs(miss.dL_abs - 180.0 / 67.0) < 1e-9);
        check("missing_spin_cross_rel", fabs(miss.rel_gross - 60.0 / 67.0) < 1e-9);
    }

    /* ---- lever-arm precision reverse solve, checked at the boundary ---- */
    am_report h1 = run_case("arm_error_1.0mm", AM_TIER_F64,
                            am_v(0, 0.001, 0), am_v(0, 0, 0), 1);
    check("one_mm_adds_about_2permille", fabs(h1.rel_gross - sol0.j * 0.001 / 3.0) < 1e-12);
    check("one_mm_violates", !h1.pass); /* ~1.99e-3 > 1e-3 */

    am_report h01 = run_case("arm_error_0.1mm", AM_TIER_F64,
                              am_v(0, 0.0001, 0), am_v(0, 0, 0), 1);
    check("tenth_mm_passes", h01.pass && h01.rel_gross < 2.1e-4);

    am_report h04 = run_case("arm_error_0.4mm", AM_TIER_F64,
                              am_v(0, 0.0004, 0), am_v(0, 0, 0), 0);
    am_report h06 = run_case("arm_error_0.6mm", AM_TIER_F64,
                              am_v(0, 0.0006, 0), am_v(0, 0, 0), 1);
    check("point4mm_inside", h04.pass);
    check("point6mm_outside", !h06.pass);

    /* ---- float32 full arm must already meet the one-per-mille gate ---- */
    am_report f32 = run_case("full_arm_f32", AM_TIER_F32,
                             am_v(0, 0, 0), am_v(0, 0, 0), 1);
    check("f32_meets_1e3", f32.pass);

    /* ---- broad randomized eccentric battery, exact double ---- */
    int n = 20000, tested = 0, near_cancel = 0;
    double worst = 0.0, worst_total = 0.0;
    int k;
    for (k = 0; k < n; ++k) {
        am_body X, Y; am_contact cc;
        memset(&X, 0, sizeof(X)); memset(&Y, 0, sizeof(Y));
        cc.p = am_v(0, 0, 0);
        /* random unit normal */
        double nx = frand(-1, 1), ny = frand(-1, 1), nz = frand(-1, 1);
        double nl = sqrt(nx * nx + ny * ny + nz * nz) + 1e-12;
        cc.n = am_v(nx / nl, ny / nl, nz / nl);
        cc.e = frand(0.0, 1.0);
        /* random eccentric arms, place COMs at p - r */
        am_vec3 rA = am_v(frand(-0.6, 0.6), frand(-0.6, 0.6), frand(-0.6, 0.6));
        am_vec3 rB = am_v(frand(-0.6, 0.6), frand(-0.6, 0.6), frand(-0.6, 0.6));
        X.x = am_scale(rA, -1.0);
        Y.x = am_scale(rB, -1.0);
        X.m = frand(0.5, 5.0); Y.m = frand(0.5, 5.0);
        X.I = frand(0.02, 0.6); Y.I = frand(0.02, 0.6);
        X.w = am_v(frand(-6, 6), frand(-6, 6), frand(-6, 6));
        Y.w = am_v(frand(-6, 6), frand(-6, 6), frand(-6, 6));
        X.v = am_v(frand(-6, 6), frand(-6, 6), frand(-6, 6));
        Y.v = am_v(frand(-6, 6), frand(-6, 6), frand(-6, 6));
        /* force an approaching normal contact speed */
        double vn = am_dot(am_sub(Y.v, X.v), cc.n);
        double want = -frand(1.0, 14.0);
        Y.v = am_add(Y.v, am_scale(cc.n, want - vn));

        am_snapshot sp = am_snapshot_make(&X, &Y);
        if (sp.gross < 1e-3) continue;
        am_solution s = am_solve_apply(&X, &Y, &cc, AM_TIER_F64,
                                       am_v(0, 0, 0), am_v(0, 0, 0));
        if (s.status != 1) continue;
        am_report rr = am_check(&sp, &X, &Y, s.eccA, s.eccB, AM_TOL_ONE_PER_MILLE);
        tested++;
        if (rr.rel_gross > worst) worst = rr.rel_gross;
        if (rr.near_cancel) near_cancel++;
        else if (rr.rel_total > worst_total) worst_total = rr.rel_total;
        if (!rr.pass) emit("random_unexpected", &rr);
    }
    printf("random_battery_cases=%d worst_rel=%.3e worst_literal_total=%.3e near_cancel=%d\n",
           tested, worst, worst_total, near_cancel);
    check("random_battery_nonempty", tested > 1000);
    check("random_all_within_1e3", worst < 1e-9);
    /* literal |dL|/|L_before| must also hold whenever the denominator is
     * meaningful (not physically cancelled toward zero). */
    check("random_literal_total_ratio", worst_total < 1e-8);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
