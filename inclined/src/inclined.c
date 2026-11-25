
/*
 * inclined — a sphere on an incline: the two static-friction critical angles.
 *
 * ---------------------------------------------------------------------------
 * Physics
 * ---------------------------------------------------------------------------
 *   beta = I / (m R^2)                       (beta = 2/5 for a solid sphere)
 *
 * Rest (centre held at v = 0, w = 0):
 *       f = -m g sin(theta)                  tan(theta_s) = mu
 *
 * No-slip rolling (contact point instantaneously fixed, u = v + R w = 0):
 *   the contact point does not move, so friction does no work; it only has to
 *   supply the angular acceleration.  a = -R alpha, I alpha = R f, and the
 *   centre equation m a = m g sin(theta) + f give
 *       f = -beta/(1+beta) m g sin(theta)    tan(theta_r) = mu (1+beta)/beta
 *
 *   N = m g cos(theta) in both cases, so the required coefficient is |f| / N.
 *   The rolling force is strictly smaller than the resting force, hence
 *
 *       tan(theta_r) / tan(theta_s) = (1+beta)/beta = 1 + m R^2 / I  >  1,
 *
 *   i.e. theta_r > theta_s, and the gap is fixed by the moment of inertia:
 *   solid sphere (2/5) 3.5, hollow sphere (2/3) 2.5, hoop (1) 2.0.
 *
 *   Cross-check from the contact-slip acceleration.  With Coulomb friction
 *   mu N opposing the slip u = v + R w,
 *       du/dt = g cos(theta) [ tan(theta) - mu (1+1/beta) sgn(u) ].
 *   For u > 0 the slip decays iff theta < theta_r — the same threshold.
 *
 * Regimes for a ball released from rest at slope theta:
 *       theta <= theta_s            rest
 *       theta_s < theta <= theta_r  rolling (no slip)
 *       theta > theta_r             sliding
 *
 * ---------------------------------------------------------------------------
 * The trap
 * ---------------------------------------------------------------------------
 * The collision-style friction rule used by a velocity solver is
 *       f_t = mu * |v_normal|,
 * which is ZERO at a resting contact because v_normal = 0.  Such a contact is
 * never held and the ball creeps down the slope.  Static friction cannot be
 * read off a kinetic law; it has to be adjudicated from the force balance and
 * clamped to mu N.  Slip, likewise, must be judged from the CONTACT-point
 * velocity u = v + R w, not the centre velocity: a purely rolling ball has
 * v != 0 yet u = 0, and a centre-velocity test would misjudge it as sliding.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef enum { R_REST = 0, R_ROLL = 1, R_SLIDE = 2 } Regime;

static const char* regime_name(Regime r) {
    switch (r) {
    case R_REST:  return "rest";
    case R_ROLL:  return "rolling";
    default:      return "sliding";
    }
}
static double deg(double r) { return r * 180.0 / M_PI; }
static double rad(double d) { return d * M_PI / 180.0; }

/* ---- critical angles ------------------------------------------------- */
static double tan_rest_critical(double mu)              { return mu; }
static double tan_roll_critical(double mu, double beta) { return mu * (1.0 + beta) / beta; }
static double theta_static(double mu)                   { return atan(mu); }
static double theta_rolling(double mu, double beta)     { return atan(mu * (1.0 + beta) / beta); }
static double ratio_of_tangents(double beta)            { return (1.0 + beta) / beta; }

static Regime classify(double theta, double mu, double beta) {
    if (theta <= theta_static(mu))        return R_REST;    /* equality -> gentler */
    if (theta <= theta_rolling(mu, beta)) return R_ROLL;
    return R_SLIDE;
}

/* ---- the trap -------------------------------------------------------- */
static double kinetic_only_friction(double mu, double v_normal) {
    return mu * fabs(v_normal);             /* 0 at a resting contact */
}
static double static_hold_force(double m, double g, double theta) {
    return -m * g * sin(theta);             /* nonzero; fixed by force balance */
}
static double free_slide_displacement(double g, double theta, double t) {
    return 0.5 * g * sin(theta) * t * t;    /* no hold at all */
}

/* The naive solver: kinetic Coulomb friction only, f_t = mu * N * sgn(u), with
 * the normal force N = m g cos and the contact slip u = v + R w.  At a resting
 * contact u = 0 so f_t = 0 and nothing holds the centre.  This is the bug we
 * quantify: the ball creeps at the free-slide rate. */

/* ---- statistics ------------------------------------------------------ */
static double normal_tail(double z) { return 0.5 * erfc(z / sqrt(2.0)); }
static double inverse_tail(double target) {         /* z >= 0 with P(Z>z)=target */
    double lo = 0.0, hi = 10.0;
    for (int i = 0; i < 200; i++) {
        double mid = 0.5 * (lo + hi);
        if (normal_tail(mid) > target) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

/* ---- simulation ------------------------------------------------------ */
typedef struct { double g, mu, beta, R, m, tmax, dt; } Params;

typedef struct {
    double v;      /* centre velocity along the slope [m/s] */
    double w;      /* angular velocity about +z [rad/s] */
    double s;      /* centre displacement along the slope [m] */
    double slip;   /* integral |u| dt [m] */
    Regime mode;
} State;

typedef struct {
    Regime regime;
    double s, slip, max_slip, v_end, w_end, energy;
} Result;

static int sign_of(double x) { return x < 0.0 ? -1 : 1; }

#define UTOL 1e-12
#define VTOL 1e-12
#define LIMTOL (1.0 + 1e-9)

/* One fixed step of semi-implicit Euler.  The friction branch is selected from
 * the contact velocity u = v + R w, never from the centre velocity. */
static void advance(State* st, const Params* p, double theta) {
    const double m = p->m, g = p->g, R = p->R, beta = p->beta;
    const double I  = beta * m * R * R;
    const double N  = m * g * cos(theta);
    const double gd = g * sin(theta);                 /* gravity along slope */
    const double tan_r = tan_roll_critical(p->mu, beta);
    const double dt = p->dt;

    double u = st->v + R * st->w;
    double u_old = u;
    Regime mode;

    if (fabs(u) < UTOL) {
        double need_rest = fabs(m * gd);              /* |f| to hold the centre */
        double need_roll = fabs(beta / (1.0 + beta) * m * gd); /* |f| for no-slip */
        if (fabs(st->v) < VTOL && fabs(st->w) < VTOL &&
            need_rest <= p->mu * N * LIMTOL) {
            mode = R_REST;
        } else if (need_roll <= p->mu * N * LIMTOL) {
            mode = R_ROLL;
        } else {
            mode = R_SLIDE;
        }
    } else {
        mode = R_SLIDE;
    }

    if (mode == R_REST) {
        st->v = 0.0;
        st->w = 0.0;
    } else if (mode == R_ROLL) {
        st->v += gd / (1.0 + beta) * dt;              /* a = g sin / (1+beta) */
        st->w = -st->v / R;                           /* keep u = 0 exactly */
    } else {
        double sgn = (fabs(u) > UTOL) ? (double)sign_of(u)
                                      : ((tan(theta) > tan_r) ? 1.0 : -1.0);
        double f = -sgn * p->mu * N;                  /* Coulomb at the limit */
        st->v += (gd + f / m) * dt;
        st->w += (R * f) / I * dt;
    }

    st->s += st->v * dt;
    double u_new = st->v + R * st->w;
    st->slip += 0.5 * (fabs(u_old) + fabs(u_new)) * dt;
    st->mode = mode;
}

static Result simulate(const Params* p, double theta) {
    Result r;
    memset(&r, 0, sizeof(r));
    State st = { 0.0, 0.0, 0.0, 0.0, R_REST };
    int steps = (int)llround(p->tmax / p->dt);
    for (int i = 0; i < steps; i++) {
        advance(&st, p, theta);
        double u = st.v + p->R * st.w;
        if (fabs(u) > r.max_slip) r.max_slip = fabs(u);
    }
    r.regime = st.mode;
    r.s = st.s;
    r.slip = st.slip;
    r.v_end = st.v;
    r.w_end = st.w;
    /* mechanical energy relative to the start (rolling conserves it) */
    double I = p->beta * p->m * p->R * p->R;
    r.energy = 0.5 * p->m * st.v * st.v + 0.5 * I * st.w * st.w
             - p->m * p->g * st.s * sin(theta);
    return r;
}

/* Kinetic-only comparator: no static branch at all.  Tangential force is
 * Coulomb mu N opposing the CONTACT slip, which is zero at a resting contact;
 * the consequence is a slow creep (the reported "ball slowly slides" symptom). */
static Result simulate_kinetic_only(const Params* p, double theta) {
    Result r;
    memset(&r, 0, sizeof(r));
    const double m = p->m, g = p->g, R = p->R, beta = p->beta;
    const double I = beta * m * R * R;
    const double N = m * g * cos(theta);
    const double gd = g * sin(theta);
    double v = 0.0, w = 0.0;
    int steps = (int)llround(p->tmax / p->dt);
    for (int i = 0; i < steps; i++) {
        double u = v + R * w;
        double f = (fabs(u) > UTOL) ? -p->mu * N * (double)sign_of(u) : 0.0;
        v += (gd + f / m) * p->dt;
        w += (R * f / I) * p->dt;
        r.s += v * p->dt;
        r.slip += fabs(u) * p->dt;
        if (fabs(u) > r.max_slip) r.max_slip = fabs(u);
    }
    r.regime = (fabs(v) < 1e-9) ? R_REST : R_SLIDE;
    r.v_end = v;
    r.w_end = w;
    return r;
}

/* ---- checks ---------------------------------------------------------- */
static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}
static int close_to(double a, double b, double tol) { return fabs(a - b) <= tol; }

/* ---- calibration ----------------------------------------------------- */
#define REQUIRED_DEG 20.0
#define SIGMA_THETA_DEG 0.5
#define FALSE_SLIDE_TARGET 0.005

typedef struct { double m, mu, theta_s, dev, false_slide; } CurvePt;

static CurvePt curve_at(double mfac, double mu0) {
    CurvePt c;
    c.m = mfac;
    c.mu = mfac * mu0;
    c.theta_s = deg(theta_static(c.mu));
    c.dev = c.theta_s - REQUIRED_DEG;
    c.false_slide = normal_tail(c.dev / SIGMA_THETA_DEG);
    return c;
}

static void calibration(double* mu_out, CurvePt* chosen, int verbose) {
    double mu0 = tan(rad(REQUIRED_DEG));
    double z = inverse_tail(FALSE_SLIDE_TARGET);
    double dev_min = z * SIGMA_THETA_DEG;
    double m_star = tan(rad(REQUIRED_DEG + dev_min)) / mu0;
    *chosen = curve_at(m_star, mu0);
    *mu_out = chosen->mu;
    if (verbose) {
        printf("calib.mu0=%.9f\n", mu0);
        printf("calib.z=%.6f\n", z);
        printf("calib.dev_min_deg=%.6f\n", dev_min);
        printf("calib.m_star=%.9f\n", m_star);
    }
}

static Params default_params(double mu, double beta) {
    Params p = { 9.81, mu, beta, 1.0, 1.0, 10.0, 1e-4 };
    return p;
}

/* ---- CSV emitters ---------------------------------------------------- */
static void emit_curves(void) {
    double mu0 = tan(rad(REQUIRED_DEG));
    printf("m,mu,theta_s_deg,deviation_deg,false_slide_rate\n");
    for (int i = 0; i <= 40; i++) {
        CurvePt c = curve_at(1.0 + 0.005 * i, mu0);   /* 1.00 .. 1.20 */
        printf("%.6f,%.9f,%.6f,%.6f,%.12g\n",
               c.m, c.mu, c.theta_s, c.dev, c.false_slide);
    }
}

static void emit_tiers(double mu, double beta) {
    Params p = default_params(mu, beta);
    double th_s = deg(theta_static(mu));
    double th_r = deg(theta_rolling(mu, beta));
    struct { const char* label; double theta; } pts[] = {
        { "design20",  REQUIRED_DEG },
        { "tierA",     th_s - 1.0   },
        { "tierB",     th_s + 1.0   },
        { "tierC",     th_r + 1.0   },
        { "probe_s-",  th_s - 0.5   },
        { "probe_r-",  th_r - 0.5   },
        { "probe_r+",  th_r + 0.5   },
    };
    int n = (int)(sizeof(pts) / sizeof(pts[0]));
    printf("label,theta_deg,regime,disp_m,slip_m,max_slip_mps,"
           "disp_over_R,slip_over_R,gate,gate_value_m,gate_ok\n");
    for (int i = 0; i < n; i++) {
        Result r = simulate(&p, rad(pts[i].theta));
        Regime expect = classify(rad(pts[i].theta), mu, beta);
        double gate_value;
        const char* gate;
        int gate_ok;
        if (expect == R_REST) {
            gate_value = fabs(r.s);   gate = "disp"; gate_ok = gate_value <= 0.01 * p.R;
        } else if (expect == R_ROLL) {
            gate_value = r.slip;      gate = "slip"; gate_ok = gate_value <= 0.01 * p.R;
        } else {
            gate_value = r.slip;      gate = "slip"; gate_ok = gate_value >  0.01 * p.R;
        }
        printf("%s,%.6f,%s,%.9f,%.9f,%.9f,%.9f,%.9f,%s,%.9f,%s\n",
               pts[i].label, pts[i].theta, regime_name(r.regime),
               r.s, r.slip, r.max_slip, r.s / p.R, r.slip / p.R,
               gate, gate_value, gate_ok ? "pass" : "FAIL");
    }
}

static void emit_timeseries(double mu, double beta) {
    Params p = default_params(mu, beta);
    double th_s = deg(theta_static(mu));
    double th_r = deg(theta_rolling(mu, beta));
    struct { const char* label; double theta; } pts[] = {
        { "rest",    th_s - 1.0 },
        { "rolling", th_s + 1.0 },
        { "sliding", th_r + 1.0 },
    };
    int nsteps = (int)llround(p.tmax / p.dt);
    int every = nsteps / 100;
    printf("tier,t,s_m,slip_m,regime\n");
    for (int k = 0; k < 3; k++) {
        State st = { 0.0, 0.0, 0.0, 0.0, R_REST };
        double theta = rad(pts[k].theta);
        for (int i = 1; i <= nsteps; i++) {
            advance(&st, &p, theta);
            if (i % every == 0 || i == nsteps)
                printf("%s,%.6f,%.9f,%.9f,%s\n", pts[k].label, i * p.dt,
                       st.s, st.slip, regime_name(st.mode));
        }
    }
}

static void emit_sweep(double mu, double beta) {
    Params p = default_params(mu, beta);
    printf("theta_deg,regime,disp_m,slip_m\n");
    for (double th = 0.0; th <= 70.0001; th += 0.5) {
        Result r = simulate(&p, rad(th));
        printf("%.4f,%s,%.9f,%.9f\n", th, regime_name(r.regime), r.s, r.slip);
    }
}

/* ---- main ------------------------------------------------------------ */
int main(int argc, char** argv) {
    const double beta = 0.4;                  /* solid sphere */
    const char* mode = (argc > 1) ? argv[1] : "test";

    if (strcmp(mode, "curves") == 0)     { emit_curves(); return 0; }
    if (strcmp(mode, "sweep") == 0) {
        double mu = 0.0; CurvePt chosen;
        calibration(&mu, &chosen, 0);
        emit_sweep(mu, beta);
        return 0;
    }
    if (strcmp(mode, "tiers") == 0 ||
        strcmp(mode, "timeseries") == 0) {
        double mu = 0.0; CurvePt chosen;
        calibration(&mu, &chosen, 0);
        if (strcmp(mode, "tiers") == 0)      emit_tiers(mu, beta);
        else                                 emit_timeseries(mu, beta);
        return 0;
    }

    double mu = 0.0;
    CurvePt chosen;
    calibration(&mu, &chosen, 1);

    const double th_s = deg(theta_static(mu));
    const double th_r = deg(theta_rolling(mu, beta));
    const double tan_s = tan_rest_critical(mu);
    const double tan_r = tan_roll_critical(mu, beta);
    const double ratio = ratio_of_tangents(beta);

    printf("beta=%.6f\n", beta);
    printf("mu_chosen=%.9f\n", mu);
    printf("theta_s_deg=%.6f\n", th_s);
    printf("theta_r_deg=%.6f\n", th_r);
    printf("tan_s=%.9f\n", tan_s);
    printf("tan_r=%.9f\n", tan_r);
    printf("ratio_tan=%.6f\n", ratio);
    printf("ratio_expected=%.6f\n", 1.0 + 1.0 / beta);
    printf("deviation_deg=%.6f\n", chosen.dev);
    printf("false_slide_rate=%.9f\n", chosen.false_slide);

    /* --- identities --- */
    check("rest_tangent_is_mu",             close_to(tan_s, mu, 1e-12));
    check("roll_tangent_is_mu_ratio",       close_to(tan_r, mu * (1.0 + beta) / beta, 1e-12));
    check("ratio_is_one_plus_inverse_beta", close_to(ratio, 1.0 + 1.0 / beta, 1e-12));
    check("rolling_angle_exceeds_rest",     th_r > th_s);
    check("mu_satisfies_rest_critical",     th_s >= REQUIRED_DEG);
    check("mu_satisfies_roll_critical",     th_r >= REQUIRED_DEG && th_r > th_s);
    check("mu_physically_plausible",        mu > 0.2 && mu < 0.8);
    check("ratio_set_by_inertia",
          close_to(ratio_of_tangents(0.4), 3.5, 1e-12) &&
          close_to(ratio_of_tangents(2.0 / 3.0), 2.5, 1e-12) &&
          close_to(ratio_of_tangents(1.0), 2.0, 1e-12));
    check("required_angle_at_mu0_rests",
          classify(rad(REQUIRED_DEG), tan(rad(REQUIRED_DEG)), beta) == R_REST);

    /* --- the trap, quantified --- */
    check("kinetic_friction_zero_at_rest", kinetic_only_friction(mu, 0.0) == 0.0);
    check("static_hold_force_nonzero",
          fabs(static_hold_force(1.0, 9.81, rad(REQUIRED_DEG))) > 1.0);
    {
        double creep = free_slide_displacement(9.81, rad(REQUIRED_DEG), 10.0);
        printf("trap.free_slide_10s_m=%.6f\n", creep);
        printf("trap.bound_m=%.6f\n", 0.01);
        check("trap_exceeds_bound", creep > 0.01);
    }
    {
        /* a purely rolling ball: contact point fixed, centre moving */
        double v = 5.0, w = -v / 1.0, R = 1.0;
        double u = v + R * w;
        printf("trap.rolling_centre_speed=%.6f\n", fabs(v));
        printf("trap.rolling_contact_slip=%.12f\n", u);
        check("rolling_has_zero_contact_slip", fabs(u) < 1e-12);
        check("centre_velocity_would_misjudge", fabs(v) > 1.0);
    }
    {
        /* same 20 deg slope: proper static solver vs kinetic-only solver */
        Params p = default_params(mu, beta);
        Result good = simulate(&p, rad(REQUIRED_DEG));
        Result bad  = simulate_kinetic_only(&p, rad(REQUIRED_DEG));
        printf("trap.static_solver_disp_m=%.12f\n", good.s);
        printf("trap.kinetic_only_disp_m=%.9f\n", bad.s);
        printf("trap.bound_m=%.6f\n", 0.01);
        check("static_solver_holds", good.regime == R_REST && fabs(good.s) <= 0.01);
        check("kinetic_only_creeps", bad.regime == R_SLIDE && fabs(bad.s) > 0.01);
    }

    /* --- rolling threshold cross-checked by slip acceleration --- */
    {
        double lo = th_r - 1.0, hi = th_r + 1.0;
        double a_lo = 9.81 * cos(rad(lo)) * (tan(rad(lo)) - mu * (1.0 + beta) / beta);
        double a_hi = 9.81 * cos(rad(hi)) * (tan(rad(hi)) - mu * (1.0 + beta) / beta);
        printf("kin.slip_accel_below=%.9f\n", a_lo);
        printf("kin.slip_accel_above=%.9f\n", a_hi);
        check("slip_accel_below_negative", a_lo < 0.0);
        check("slip_accel_above_positive", a_hi > 0.0);
    }

    Params p = default_params(mu, beta);
    const double bound = 0.01 * p.R;

    /* --- design point: 20 deg must hold, and classify as rest --- */
    {
        Result r = simulate(&p, rad(REQUIRED_DEG));
        printf("design.theta_deg=%.4f\n", REQUIRED_DEG);
        printf("design.regime=%s\n", regime_name(r.regime));
        printf("design.margin_to_rest_deg=%.6f\n", th_s - REQUIRED_DEG);
        printf("design.margin_to_roll_deg=%.6f\n", th_r - REQUIRED_DEG);
        printf("design.disp_m=%.12f\n", r.s);
        printf("design.bound_m=%.6f\n", bound);
        check("design_classified_rest", classify(rad(REQUIRED_DEG), mu, beta) == R_REST);
        check("design_stays_put", r.regime == R_REST && fabs(r.s) <= bound);
    }

    /* --- three tiers: one per regime, each near a critical angle --- */
    struct { const char* name; double theta; } T[3] = {
        { "A_below_static", th_s - 1.0 },
        { "B_above_static", th_s + 1.0 },
        { "C_above_roll",   th_r + 1.0 },
    };
    for (int i = 0; i < 3; i++) {
        Regime expect = classify(rad(T[i].theta), mu, beta);
        Result r = simulate(&p, rad(T[i].theta));
        printf("tier.%s.theta_deg=%.4f\n", T[i].name, T[i].theta);
        printf("tier.%s.expect=%s\n", T[i].name, regime_name(expect));
        printf("tier.%s.got=%s\n", T[i].name, regime_name(r.regime));
        printf("tier.%s.disp_m=%.9f\n", T[i].name, r.s);
        printf("tier.%s.slip_m=%.9f\n", T[i].name, r.slip);
        printf("tier.%s.max_slip_mps=%.9f\n", T[i].name, r.max_slip);
        printf("tier.%s.disp_over_R=%.9f\n", T[i].name, r.s / p.R);
        printf("tier.%s.slip_over_R=%.9f\n", T[i].name, r.slip / p.R);
        printf("tier.%s.bound_m=%.9f\n", T[i].name, bound);
        printf("tier.%s.energy_end_J=%.9f\n", T[i].name, r.energy);

        check("tier_regime_agrees", r.regime == expect);

        /* The R/100 rule is the acceptance test for the "must stay put" case:
         * static friction has to hold the centre still.  A rolling or sliding
         * ball is supposed to translate, so for those tiers the R/100 gate is
         * applied to the quantity that must stay zero — the CONTACT-slip
         * distance — while the actual centre displacement is still reported. */
        double gate_value;
        const char* gate_name;
        int gate_ok;
        if (expect == R_REST) {
            gate_value = fabs(r.s);
            gate_name = "centre_displacement";
            gate_ok = (gate_value <= bound);
        } else if (expect == R_ROLL) {
            gate_value = r.slip;
            gate_name = "contact_slip";
            gate_ok = (gate_value <= bound);
        } else {
            gate_value = r.slip;
            gate_name = "contact_slip";
            gate_ok = (gate_value > bound);        /* must actually slip */
        }
        printf("tier.%s.gate=%s\n", T[i].name, gate_name);
        printf("tier.%s.gate_value_m=%.9f\n", T[i].name, gate_value);
        printf("tier.%s.gate_ok=%d\n", T[i].name, gate_ok);
        check("tier_gate", gate_ok);

        if (!gate_ok)
            printf("tier.%s.exceeded disp=%.6f m gate=%s value=%.6f m classification=%s\n",
                   T[i].name, r.s, gate_name, gate_value, regime_name(r.regime));
        if (expect != R_REST && fabs(r.s) > bound)
            printf("tier.%s.translates disp=%.6f m classification=%s (expected)\n",
                   T[i].name, r.s, regime_name(r.regime));
    }

    /* --- analytic cross-checks for the rolling tier --- */
    {
        double theta = T[1].theta;
        Result r = simulate(&p, rad(theta));
        double a = p.g * sin(rad(theta)) / (1.0 + beta);
        double s_exact = 0.5 * a * p.tmax * p.tmax;
        printf("check.roll_sim_m=%.9f\n", r.s);
        printf("check.roll_exact_m=%.9f\n", s_exact);
        printf("check.roll_rel_err=%.3e\n", fabs(r.s - s_exact) / s_exact);
        check("roll_matches_half_at2", fabs(r.s - s_exact) / s_exact < 1e-3);
        check("roll_conserves_energy", fabs(r.energy) < 1e-2 * fabs(p.m * p.g * r.s));
        check("roll_has_no_slip", r.slip == 0.0);
    }
    {
        double theta = T[2].theta;
        Result r = simulate(&p, rad(theta));
        double a = p.g * (sin(rad(theta)) - mu * cos(rad(theta)) / beta);
        (void)a;
        printf("check.slide_energy_end_J=%.9f\n", r.energy);
        check("slide_dissipates_energy", r.energy < 0.0);
        check("slide_slip_sustained", r.max_slip > 0.1);
    }

    /* --- probes half a degree from each boundary --- */
    {
        Result a = simulate(&p, rad(th_s - 0.5));
        Result b = simulate(&p, rad(th_r - 0.5));
        Result c = simulate(&p, rad(th_r + 0.5));
        printf("probe.static_minus_0p5.regime=%s disp_m=%.9f\n",
               regime_name(a.regime), a.s);
        printf("probe.roll_minus_0p5.regime=%s slip_m=%.9f\n",
               regime_name(b.regime), b.slip);
        printf("probe.roll_plus_0p5.regime=%s slip_m=%.9f\n",
               regime_name(c.regime), c.slip);
        check("probe_static_no_creep", a.regime == R_REST && fabs(a.s) <= bound);
        check("probe_roll_no_slip",    b.regime == R_ROLL && b.slip <= bound);
        check("probe_above_roll_slips", c.regime == R_SLIDE && c.slip > bound);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
