
/* order_compare: normal-first vs tangent-first friction budget ordering.
 *
 * Unit-mass sphere vs plane, n orthogonal to t:
 *   v^- = -u*n + s*t, u > 0, s >= 0.
 *
 * The program emits CSV for the five fixed oblique impacts used in
 * artifacts/friction_order/report.md.  The Python driver owns the sweep and
 * plot; this binary is the independent C cross-check of the five rows.
 */
#include <stdio.h>
#include <math.h>

#define RUN_MU 0.30

typedef struct {
    int id;
    double u, s, e, mu;
} Case;

static double clamp_min(double x, double lo) { return x < lo ? lo : x; }
static double min_d(double a, double b) { return a < b ? a : b; }

int main(void) {
    const Case cases[5] = {
        {1, 10.0, 8.0, 0.50, RUN_MU},
        {2, 10.0, 7.0, 0.60, RUN_MU},
        {3, 10.0, 6.5, 0.75, RUN_MU},
        {4, 10.0, 6.1, 0.90, RUN_MU},
        {5, 10.0, 3.0, 1.00, RUN_MU},
    };

    puts("id,order,u,s,e,mu,dt,vt_after,jn_per_mass,cone_bound,utilization,excess_pct,cone_ok,projected_dt,projected_vt,projected_ok");

    for (int i = 0; i < 5; ++i) {
        const Case *c = &cases[i];
        double jn = (1.0 + c->e) * c->u;
        double bound = c->mu * jn;

        /* Normal-first: friction uses the frozen uncorrected incoming load u. */
        double n_dt = min_d(c->s, c->mu * c->u);
        double n_use = n_dt / bound;
        double n_excess = 100.0 * clamp_min(n_use - 1.0, 0.0);
        int n_ok = n_dt <= bound + 1e-10;
        double n_proj = min_d(n_dt, bound);

        /* Tangent-first: before jn exists, use the maximal reversal 2u.
         * For e < 1 this can be wider than the true Coulomb cone. */
        double t_dt = min_d(c->s, 2.0 * c->mu * c->u);
        double t_use = t_dt / bound;
        double t_excess = 100.0 * clamp_min(t_use - 1.0, 0.0);
        int t_ok = t_dt <= bound + 1e-10;
        double t_proj = min_d(t_dt, bound);

        printf("%d,normal_first,%.17g,%.17g,%.17g,%.17g,%.6f,%.6f,%.6f,%.6f,%.9f,%.6f,%d,%.6f,%.6f,%d\n",
               c->id, c->u, c->s, c->e, c->mu,
               n_dt, c->s - n_dt, jn, bound, n_use, n_excess, n_ok,
               n_proj, c->s - n_proj, n_proj <= bound + 1e-10);
        printf("%d,tangent_first,%.17g,%.17g,%.17g,%.17g,%.6f,%.6f,%.6f,%.6f,%.9f,%.6f,%d,%.6f,%.6f,%d\n",
               c->id, c->u, c->s, c->e, c->mu,
               t_dt, c->s - t_dt, jn, bound, t_use, t_excess, t_ok,
               t_proj, c->s - t_proj, t_proj <= bound + 1e-10);
    }
    return 0;
}
