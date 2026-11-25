
/* chain: one-dimensional elastic collision chain.
 *
 * A row of balls, only the first moving.  With increasing masses each ball
 * recoils after its rear collision, so every neighbouring pair collides
 * exactly once and the chain is a sequence of two-body collisions.
 *
 * Derivation, one elastic link (ball k, mass m1, speed u, hits ball k+1 at
 * rest, mass m2):
 *   momentum:      m1 u        = m1 v1 + m2 v2
 *   energy:    1/2 m1 u^2      = 1/2 m1 v1^2 + 1/2 m2 v2^2
 *   => v2 = 2 m1 u / (m1 + m2)
 *   => E_next/E_before = m2 v2^2 / (m1 u^2) = 4 m1 m2 / (m1 + m2)^2.
 * The energy that reaches ball n is therefore
 *   E_n / E_1 = prod_{k=1..n-1} 4 m_k m_{k+1} / (m_k + m_{k+1})^2,
 * because ball k first receives E_k from ball k-1 and then splits it with
 * ball k+1 in the same ratio.  With a constant growth ratio
 * rho = m_{k+1}/m_k this is a geometric product:
 *   T(rho)   = 4 rho / (1 + rho)^2                (per-link efficiency)
 *   E_n/E_1  = T^(n-1)
 *   L(rho)   = ln(0.1) / ln T                     links until 10% survives
 *   rho(L)   = (2 - T + 2 sqrt(1-T)) / T,         T = 0.1^(1/L)
 * The growing branch rho >= 1 is the one used below.
 *
 * The front balls recoil with |v_1| = (rho-1)/(rho+1) * u0 and
 * |v_2| = 2(rho-1)/(rho+1)^2 * u0 = q |v_1|, q = 2/(1+rho); these are the
 * fly-off speeds that trade against a short L.
 *
 * The self-test checks the formulas, the inverse design and an event-driven
 * simulation of the chain against the analytic velocities.
 */
#include <stdio.h>
#include <math.h>

#define MAXN 8

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* ---- closed form ------------------------------------------------------ */

/* fraction of the moving ball's energy that lands in the next ball */
static double link_energy_frac(double m1, double m2) {
    double s = m1 + m2;
    return 4.0 * m1 * m2 / (s * s);
}

/* same for a growth ratio rho = m_next / m_prev */
static double ratio_energy_frac(double rho) {
    double s = 1.0 + rho;
    return 4.0 * rho / (s * s);
}

/* inverse of ratio_energy_frac on the growing branch (rho >= 1) */
static double ratio_for_frac(double t) {
    return ((2.0 - t) + 2.0 * sqrt(1.0 - t)) / t;
}

/* number of links until the surviving energy fraction reaches 10% */
static double propagation_length(double rho) {
    if (rho <= 1.0) return INFINITY;          /* equal masses: no decay */
    return log(0.1) / log(ratio_energy_frac(rho));
}

static double ratio_for_length(double L) {
    return ratio_for_frac(pow(0.1, 1.0 / L));
}

static void fill_uniform(double* m, int n, double rho) {
    m[0] = 1.0;
    for (int i = 1; i < n; ++i) m[i] = m[i - 1] * rho;
}

/* ---- analytic chain --------------------------------------------------- */

/* v: final velocities, recv: velocity the ball had right after being hit */
static void chain_analytic(const double* m, int n, double u0,
                           double* v, double* recv) {
    for (int i = 0; i < n; ++i) { v[i] = 0.0; recv[i] = 0.0; }
    recv[0] = u0;
    double u = u0;
    for (int k = 0; k + 1 < n; ++k) {
        double s = m[k] + m[k + 1];
        v[k] = (m[k] - m[k + 1]) / s * u;
        u = 2.0 * m[k] / s * u;
        v[k + 1] = u;
        recv[k + 1] = u;
    }
}

/* ---- event-driven simulation ----------------------------------------- */

/* balls of diameter 1, small initial gaps; advance to the next contact,
 * apply the elastic impulse, repeat.  Increasing masses make every pair
 * collide once, so the event order reproduces the analytic chain.  The last
 * event is the one at ball n-1, which is when it reaches its final speed
 * (it never collides again), so the chain cannot reorder.  The simulation
 * is deliberately independent of chain_analytic: with equal or decreasing
 * masses it produces the multiple-collision order instead. */
static int chain_simulate(const double* m, int n, double u0, double gap,
                          double* vout, double* recvout) {
    const double dia = 1.0;
    double x[MAXN], v[MAXN];
    for (int i = 0; i < n; ++i) {
        x[i] = i * (dia + gap);
        v[i] = (i == 0) ? u0 : 0.0;
        recvout[i] = 0.0;
    }
    recvout[0] = u0;
    int events = 0;

    for (int guard = 0; guard < 10000; ++guard) {
        double best = INFINITY;
        int p = -1;
        for (int i = 0; i + 1 < n; ++i) {
            double sep = x[i + 1] - x[i] - dia;
            double rel = v[i] - v[i + 1];
            if (rel > 1e-15) {
                double dt = sep / rel;
                if (dt < 0.0) dt = 0.0;
                if (dt < best) { best = dt; p = i; }
            }
        }
        if (p < 0) break;                     /* nothing else can collide */

        for (int i = 0; i < n; ++i) x[i] += v[i] * best;
        x[p + 1] = x[p] + dia;                /* snap to exact contact */

        double m1 = m[p], m2 = m[p + 1];
        double u1 = v[p], u2 = v[p + 1];
        v[p]     = ((m1 - m2) * u1 + 2.0 * m2 * u2) / (m1 + m2);
        v[p + 1] = ((m2 - m1) * u2 + 2.0 * m1 * u1) / (m1 + m2);
        if (recvout[p + 1] == 0.0 && v[p + 1] > 0.0) recvout[p + 1] = v[p + 1];
        events++;
    }

    for (int i = 0; i < n; ++i) vout[i] = v[i];
    return events;
}

static double kinetic(const double* m, const double* v, int n) {
    double e = 0.0;
    for (int i = 0; i < n; ++i) e += 0.5 * m[i] * v[i] * v[i];
    return e;
}

static double momentum(const double* m, const double* v, int n) {
    double p = 0.0;
    for (int i = 0; i < n; ++i) p += m[i] * v[i];
    return p;
}

/* received energy / initial energy, from the recorded post-hit velocity */
static double received_ratio(const double* m, const double* recv, int k,
                             double u0) {
    return m[k] * recv[k] * recv[k] / (m[0] * u0 * u0);
}

static double max_abs_diff(const double* a, const double* b, int n) {
    double d = 0.0;
    for (int i = 0; i < n; ++i) d = fmax(d, fabs(a[i] - b[i]));
    return d;
}

/* recoil speed of the first two balls, in units of u0 */
static double flyoff_first(double rho) {
    return (rho - 1.0) / (rho + 1.0);
}

static double flyoff_second(double rho) {
    double q = 2.0 / (1.0 + rho);
    return flyoff_first(rho) * q;
}

int main(void) {
    const double u0 = 1.0, gap = 1e-3;
    double m[MAXN], va[MAXN], vs[MAXN], ra[MAXN], rs[MAXN];

    /* -- 1. single-link energy fraction matches the collision law -------- */
    {
        double M = 2.0, m2 = 6.0, u = 3.0, s = M + m2;
        double vm = 2.0 * M / s * u;
        double frac = (0.5 * m2 * vm * vm) / (0.5 * M * u * u);
        check("link_fraction_matches_collision",
              fabs(frac - link_energy_frac(M, m2)) < 1e-12);
        check("link_fraction_rho3",
              fabs(ratio_energy_frac(3.0) - 0.75) < 1e-12);
    }

    /* -- 2. inverses ----------------------------------------------------- */
    {
        int ok = 1;
        for (double r = 1.01; r < 600.0; r *= 1.37)
            if (fabs(ratio_for_frac(ratio_energy_frac(r)) - r) > 1e-9 * r) ok = 0;
        check("inverse_ratio_roundtrip", ok);

        ok = 1;
        for (double L = 1.05; L < 60.0; L *= 1.29)
            if (fabs(propagation_length(ratio_for_length(L)) - L) > 1e-9 * L) ok = 0;
        check("inverse_length_roundtrip", ok);

        check("length_equal_mass_infinite", isinf(propagation_length(1.0)));
    }

    /* -- 3. equal masses: last ball takes everything --------------------- */
    {
        fill_uniform(m, 4, 1.0);
        chain_analytic(m, 4, u0, va, ra);
        int ev = chain_simulate(m, 4, u0, gap, vs, rs);
        check("equal_mass_front_balls_stop",
              fabs(va[0]) < 1e-12 && fabs(va[1]) < 1e-12 && fabs(va[2]) < 1e-12);
        check("equal_mass_last_carries_all", fabs(va[3] - u0) < 1e-12);
        check("equal_mass_sim_matches", max_abs_diff(va, vs, 4) < 1e-9 && ev == 3);
        check("equal_mass_e4_all", fabs(received_ratio(m, ra, 3, u0) - 1.0) < 1e-12);
    }

    /* -- 4. inverse design threshold ------------------------------------- */
    /* E4/E1 < 0.02 with three links  =>  T < 0.02^(1/3)  =>  rho > 12.657 */
    double rho_min = ratio_for_frac(pow(0.02, 1.0 / 3.0));
    double L_max = propagation_length(rho_min);
    check("rho_min_for_fourth_below_2pct", fabs(rho_min - 12.6571) < 1e-3);
    check("length_at_threshold", fabs(L_max - 1.7658) < 1e-3);
    check("threshold_energy_exact",
          fabs(pow(ratio_energy_frac(rho_min), 3) - 0.02) < 1e-12);

    /* -- 5. trade-off: propagation length vs front-ball fly-off ---------- */
    printf("# tradeoff: rho T L10 d1 d2 d1+d2   (d in units of u0*tau)\n");
    {
        double probe[] = {2.0, 4.0, 6.0, 8.0, 10.0, 12.0, rho_min,
                          13.0, 15.0, 18.0, 25.0, 40.0};
        int np = (int)(sizeof(probe) / sizeof(probe[0]));
        for (int i = 0; i < np; ++i) {
            double rho = probe[i];
            double d1 = flyoff_first(rho);
            double d2 = flyoff_second(rho);
            printf("tradeoff %.4f %.6f %.4f %.4f %.4f %.4f\n",
                   rho, ratio_energy_frac(rho), propagation_length(rho),
                   d1, d2, d1 + d2);
        }
    }

    /* -- 6. three mass-ratio tiers, each simulated once ------------------ */
    printf("# tiers\n");
    {
        double tiers[3] = {13.0, 15.0, 18.0};
        for (int t = 0; t < 3; ++t) {
            double rho = tiers[t];
            fill_uniform(m, 4, rho);
            chain_analytic(m, 4, u0, va, ra);
            int ev = chain_simulate(m, 4, u0, gap, vs, rs);

            double err = max_abs_diff(va, vs, 4);
            double e4 = received_ratio(m, rs, 3, u0);
            double e4_ana = received_ratio(m, ra, 3, u0);
            double d1 = flyoff_first(rho);
            double d2 = flyoff_second(rho);

            printf("tier rho=%.3f L=%.4f d1=%.4f d2=%.4f "
                   "v=[%.6f %.6f %.6f %.6f] e4/e1=%.6f err=%.2e events=%d\n",
                   rho, propagation_length(rho), d1, d2,
                   vs[0], vs[1], vs[2], vs[3], e4, err, ev);

            char name[80];
            snprintf(name, sizeof name, "tier_%.0f_velocity_within_2pct", rho);
            check(name, err < 0.02 * u0);
            snprintf(name, sizeof name, "tier_%.0f_e4_below_2pct", rho);
            check(name, e4 < 0.02);
            snprintf(name, sizeof name, "tier_%.0f_initial_energy_conserved", rho);
            check(name, fabs(kinetic(m, vs, 4) - 0.5 * m[0] * u0 * u0) < 1e-9);
            snprintf(name, sizeof name, "tier_%.0f_momentum_conserved", rho);
            check(name, fabs(momentum(m, vs, 4) - m[0] * u0) < 1e-9);
            snprintf(name, sizeof name, "tier_%.0f_received_matches_T3", rho);
            check(name, fabs(received_ratio(m, rs, 3, u0) -
                             pow(ratio_energy_frac(rho), 3)) < 1e-9);
            snprintf(name, sizeof name, "tier_%.0f_events", rho);
            check(name, ev == 3);
            (void)e4_ana;
        }
    }

    /* -- 7. non-uniform design solved from the formula ------------------- */
    /* targets: 50% into ball 2, 10% into ball 3, 1.9% into ball 4 */
    printf("# design\n");
    {
        double target[3] = {0.50, 0.10, 0.019};
        m[0] = 1.0;
        double prev = 1.0;
        for (int k = 0; k < 3; ++k) {
            double t = target[k] / prev;       /* per-link transfer */
            m[k + 1] = m[k] * ratio_for_frac(t);
            prev = target[k];
        }
        chain_analytic(m, 4, u0, va, ra);
        int ev = chain_simulate(m, 4, u0, gap, vs, rs);

        double err = max_abs_diff(va, vs, 4);
        double e2 = received_ratio(m, rs, 1, u0);
        double e3 = received_ratio(m, rs, 2, u0);
        double e4 = received_ratio(m, rs, 3, u0);
        double prod = 1.0;
        for (int k = 0; k < 3; ++k) prod *= ratio_energy_frac(m[k + 1] / m[k]);

        printf("design masses=[%.4f %.4f %.4f %.4f] rho=[%.4f %.4f %.4f] "
               "v=[%.6f %.6f %.6f %.6f] e=[%.6f %.6f %.6f] err=%.2e events=%d\n",
               m[0], m[1], m[2], m[3],
               m[1] / m[0], m[2] / m[1], m[3] / m[2],
               vs[0], vs[1], vs[2], vs[3], e2, e3, e4, err, ev);

        check("design_mass2_solved", fabs(m[1] - 5.828427) < 1e-5);
        check("design_mass3_solved", fabs(m[2] - 104.5869) < 1e-3);
        check("design_mass4_solved", fabs(m[3] - 1987.152) < 1e-1);
        check("design_product_expr", fabs(prod - 0.019) < 1e-12);
        check("design_velocity_within_2pct", err < 0.02 * u0);
        check("design_e2_half", fabs(e2 - 0.50) < 1e-9);
        check("design_e3_tenth", fabs(e3 - 0.10) < 1e-9);
        check("design_e4_below_2pct", e4 < 0.02);
        check("design_fourth_almost_still", fabs(vs[3]) < 0.01 * u0);
        check("design_initial_energy_conserved",
              fabs(kinetic(m, vs, 4) - 0.5 * m[0] * u0 * u0) < 1e-9);
        check("design_momentum_conserved",
              fabs(momentum(m, vs, 4) - m[0] * u0) < 1e-9);
        check("design_events", ev == 3);
    }

    /* -- 8. product expression vs the direct chain for arbitrary masses -- */
    {
        double mm[4] = {1.0, 2.5, 7.0, 30.0};
        chain_analytic(mm, 4, u0, va, ra);
        double p = 1.0;
        for (int k = 0; k < 3; ++k) p *= link_energy_frac(mm[k], mm[k + 1]);
        check("product_expr_matches_chain",
              fabs(received_ratio(mm, ra, 3, u0) - p) < 1e-12);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
