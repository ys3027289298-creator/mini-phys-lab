
/* roll: rolling-sphere scenarios for the tangential-impulse investigation.
 *
 * Two reports:
 *   A) a ball rolling ONTO A SLOPE accelerates abnormally;
 *   B) a ball rolling on FLAT ground slowly veers sideways.
 *
 * The tree originally only had a point-mass resolver (solve/src/solve.c),
 * which has no angular state, so the rolling contact layer did not exist as a
 * runnable scenario.  This file adds the two reported scenarios and a small
 * residual study:
 *
 *   bad      : reported faulty tangential path -- the impulse is applied
 *              ALONG the slip direction (missing opposition sign in the
 *              shared direction calculation) and only linear velocity is
 *              updated ("rotational friction not done").
 *   dironly  : ONLY the direction/opposition sign is corrected; still no
 *              angular impulse.  This is the "can one direction fix cover
 *              both reports?" counterfactual.
 *   full     : correct contact-point slip direction
 *              u = v + omega x r, rolling effective mass
 *              k_t = 1/m + R^2/I, angular impulse domega = I^-1 (r x J),
 *              plus torsional spin friction about the contact normal.
 *
 * `make test` runs both scenarios and enforces the acceptance thresholds on
 * the full solver:
 *   A) |a_meas - 5/7 g sin theta| / a_ref <= 3%   (solid sphere, I=2/5 mR^2)
 *   B) max |z| / R                                <= 5%.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

typedef struct { double x, y, z; } V3;

static V3 vadd(V3 a, V3 b) { return (V3){a.x+b.x, a.y+b.y, a.z+b.z}; }
static V3 vsub(V3 a, V3 b) { return (V3){a.x-b.x, a.y-b.y, a.z-b.z}; }
static V3 vmul(V3 a, double s) { return (V3){a.x*s, a.y*s, a.z*s}; }
static double vdot(V3 a, V3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static V3 vcross(V3 a, V3 b) {
    return (V3){ a.y*b.z - a.z*b.y,
                 a.z*b.x - a.x*b.z,
                 a.x*b.y - a.y*b.x };
}
static double vnorm(V3 a) { return sqrt(vdot(a, a)); }

typedef struct {
    V3 p, v, w;          /* centre position, linear and angular velocity */
    double R, m, I;      /* radius, mass, scalar inertia (solid sphere) */
    double mu;           /* Coulomb friction coefficient */
    double mu_t;         /* torsional spin-friction coefficient */
} Ball;

typedef struct {
    int oppose;     /* -1: J opposes slip (correct); +1: J along slip (bug) */
    int angular;    /* apply r x J angular impulse / use rolling k_t */
    int torsion;    /* damp spin about the contact normal */
} SolverOpts;

static const SolverOpts OPT_BAD     = { +1, 0, 0 };
static const SolverOpts OPT_DIRONLY = { -1, 0, 0 };
static const SolverOpts OPT_ROTBAD  = { +1, 1, 0 };
static const SolverOpts OPT_FULL    = { -1, 1, 1 };

/* One plane contact at signed distance p.n = R, unit outward normal n. */
static void contact_step(Ball* b, V3 n, double dt, V3 gravity, SolverOpts o) {
    (void)dt;
    b->v = vadd(b->v, vmul(gravity, dt));
    b->p = vadd(b->p, vmul(b->v, dt));

    double dist = vdot(b->p, n);     /* plane through origin */
    if (dist >= b->R) return;

    V3 r = vmul(n, -b->R);           /* centre -> contact point */

    /* normal impulse, restitution 0; k_n = 1/m because r x n = 0 */
    V3 u = vadd(b->v, vcross(b->w, r));
    double un = vdot(u, n);
    double jn = 0.0;
    if (un < 0.0) {
        jn = -un * b->m;
        b->v = vadd(b->v, vmul(n, jn / b->m));
    }
    /* positional correction */
    if (vdot(b->p, n) < b->R)
        b->p = vadd(b->p, vmul(n, b->R - vdot(b->p, n)));

    /* tangential contact slip: u_t = (v + omega x r)_t */
    u = vadd(b->v, vcross(b->w, r));
    double udotn = vdot(u, n);
    V3 ut = vsub(u, vmul(n, udotn));
    double slip = vnorm(ut);

    if (slip > 1e-12 && jn > 0.0) {
        V3 shat = vmul(ut, 1.0 / slip);
        /* effective inverse mass along a tangent direction */
        double ke = 1.0 / b->m;
        if (o.angular) ke += b->R * b->R / b->I;   /* = 7/(2m) solid sphere */
        double jmag = slip / ke;                   /* impulse that stops slip */
        double cap = b->mu * jn;                   /* Coulomb cone */
        if (jmag > cap) jmag = cap;
        V3 J = vmul(shat, (double)o.oppose * jmag);
        b->v = vadd(b->v, vmul(J, 1.0 / b->m));
        if (o.angular)
            b->w = vadd(b->w, vmul(vcross(r, J), 1.0 / b->I));
    }

    /* torsional friction: kill spin omega.n, bounded by mu_t R jn */
    if (o.torsion && o.angular && jn > 0.0) {
        double wn = vdot(b->w, n);
        double Hmax = b->mu_t * b->R * jn;
        double H = fmin(fabs(b->I * wn), Hmax);
        if (fabs(wn) > 1e-12)
            b->w = vsub(b->w, vmul(n, copysign(H / b->I, wn)));
    }
}

/* least-squares slope of y vs x */
static double lin_slope(const double* x, const double* y, int n) {
    double sx=0, sy=0, sxx=0, sxy=0;
    for (int i=0;i<n;i++){ sx+=x[i]; sy+=y[i]; sxx+=x[i]*x[i]; sxy+=x[i]*y[i]; }
    double den = n*sxx - sx*sx;
    return (n*sxy - sx*sy) / den;
}

#define DT  (1.0/1000.0)
#define G   9.81
#define R_R 0.5
#define PI  3.14159265358979323846

/* Scenario A: solid sphere already in pure rolling at entry speed, then
 * running down a theta slope.  Returns fitted down-slope acceleration. */
static double slope_run(SolverOpts o, double* err_pct) {
    double theta = 20.0 * PI / 180.0;
    V3 d = { cos(theta), -sin(theta), 0.0 };   /* down-slope */
    V3 n = { sin(theta),  cos(theta), 0.0 };   /* outward normal */
    V3 l = vcross(d, n);                         /* +z */
    (void)l;
    double V0 = 2.0, T = 2.5;

    Ball b; memset(&b, 0, sizeof b);
    b.R=R_R; b.m=1.0; b.I=0.4*b.m*b.R*b.R; b.mu=0.45; b.mu_t=0.05;
    b.p = (V3){0.0, b.R, 0.0};
    b.v = vmul(d, V0);
    /* omega = -(V0/R) (d x n) gives omega x r = -V0 d (pure rolling) */
    b.w = vmul(vcross(d,n), -V0/b.R);

    V3 g = {0.0,-G,0.0};
    int N=(int)(T/DT);
    double ts[4096], vd[4096]; int k=0;
    for (int i=0;i<N;i++) {
        contact_step(&b, n, DT, g, o);
        if (i >= (int)(0.50*T/DT) && i <= (int)(0.95*T/DT)) {
            ts[k]=i*DT; vd[k]=vdot(b.v,d); k++;
        }
    }
    double a = lin_slope(ts, vd, k);
    double aref = (5.0/7.0)*G*sin(theta);       /* I=2/5 mR^2 */
    if (err_pct) *err_pct = 100.0*fabs(a-aref)/aref;
    return a;
}

/* Scenario B: straight rolling on flat ground plus a small wobble component
 * omega_x (the perturbation reports describe).  Returns max |z|/R. */
static double flat_run(SolverOpts o, double* zmax_out) {
    V3 n = {0,1,0};
    double V0=2.0, wx=0.020, T=4.0;
    Ball b; memset(&b, 0, sizeof b);
    b.R=R_R; b.m=1.0; b.I=0.4*b.m*b.R*b.R; b.mu=0.45; b.mu_t=0.05;
    b.p=(V3){0,b.R,0};
    b.v=(V3){V0,0,0};
    b.w=(V3){wx, 0.0, -V0/b.R};
    V3 g={0,-G,0};
    int N=(int)(T/DT);
    double zmax=0;
    for (int i=0;i<N;i++) {
        contact_step(&b,n,DT,g,o);
        if (fabs(b.p.z)>zmax) zmax=fabs(b.p.z);
    }
    if (zmax_out) *zmax_out=zmax;
    return 100.0*zmax/b.R;
}

/* Torsion check: a pure normal spin must decay under torsional friction. */
static double torsion_run(void) {
    V3 n={0,1,0};
    Ball b; memset(&b,0,sizeof b);
    b.R=R_R; b.m=1.0; b.I=0.4*b.m*b.R*b.R; b.mu=0.45; b.mu_t=0.2;
    b.p=(V3){0,b.R,0}; b.w=(V3){0,6.0,0};
    V3 g={0,-G,0};
    for (int i=0;i<(int)(1.0/DT);i++) contact_step(&b,n,DT,g,OPT_FULL);
    return b.w.y;
}

static int failures=0;
static void check(const char* name, int ok){
    printf("check %-28s %s\n", name, ok?"PASS":"FAIL");
    if(!ok) failures++;
}

int main(void){
    printf("mini-phys-lab roll scenarios (dt=%.0f us, solid sphere R=%.2f)\n",
           DT*1e6, R_R);
    printf("reference slope acceleration a*=5/7 g sin20 = %.4f m/s^2\n\n",
           (5.0/7.0)*G*sin(20.0*PI/180.0));

    struct { const char* tag; SolverOpts o; } modes[] = {
        {"bad_along_noRot",  OPT_BAD},
        {"dironly_noRot",    OPT_DIRONLY},
        {"wrongDir_rot",     OPT_ROTBAD},
        {"full_dirPlusRot",  OPT_FULL},
    };
    double full_aerr=0, full_zR=0;
    for (unsigned mi=0; mi<sizeof(modes)/sizeof(modes[0]); mi++){
        double aerr,zmax;
        double a=slope_run(modes[mi].o,&aerr);
        double zR=flat_run(modes[mi].o,&zmax);
        printf("[%-16s] A_slope a=%8.4f m/s2 err=%8.2f%% | "
               "B_flat zmax=%.4f m z/R=%7.3f%%\n",
               modes[mi].tag,a,aerr,zmax,zR);
        if (modes[mi].o.oppose==-1 && modes[mi].o.angular==1 &&
            modes[mi].o.torsion==1){ full_aerr=aerr; full_zR=zR; }
    }
    double wy=torsion_run();
    printf("\ntorsion: omega_y after 1 s pure spin = %.3f rad/s (initial 6)\n\n",wy);

    /* ---- acceptance thresholds on the fixed solver ---- */
    check("A_slope_acc_err_le_3pct", full_aerr <= 3.0);
    check("B_flat_lateral_le_5pctR", full_zR <= 5.0);
    check("torsional_spin_decays", fabs(wy) < 1.0);

    /* ---- lock the diagnosis: counterfactuals must really fail ---- */
    double bad_aerr, dir_aerr, dir_zR, bad_zmax;
    slope_run(OPT_BAD,&bad_aerr);
    slope_run(OPT_DIRONLY,&dir_aerr);
    dir_zR=flat_run(OPT_DIRONLY,NULL);
    (void)flat_run(OPT_BAD,&bad_zmax);
    check("diag_bad_slope_over_3pct", bad_aerr > 3.0);
    check("diag_dironly_slope_still_bad", dir_aerr > 3.0);
    check("diag_dironly_flat_still_over5", dir_zR > 5.0);

    printf("\nTHRESHOLD_A_actual_err_pct=%.3f (limit 3.000)\n", full_aerr);
    printf("THRESHOLD_B_actual_z_over_R_pct=%.3f (limit 5.000)\n", full_zR);
    printf("all_ok=%d\n", failures==0);
    return failures?1:0;
}
