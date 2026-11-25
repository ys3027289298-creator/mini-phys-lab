
/* am: collision angular-momentum budget for rigid bodies (SI units, 3D).
 *
 * Sign convention (fixed for the whole module -- do not mix any other):
 *   - one right-handed inertial world frame; +z is out of the xy plane and
 *     positive omega is counter-clockwise / right-hand-positive;
 *   - lever arm is ALWAYS r = p_contact - x_com, i.e. COM -> contact;
 *   - normal n points from body A to body B; impulse on A is -j*n and on
 *     B is +j*n (equal and opposite);
 *   - angular momentum about a single fixed origin O is
 *       L_O = sum_i [ (x_i - O) x (m_i v_i) + I_i w_i ].
 *
 * Reversing the lever-arm direction, using a screen-down y axis for omega,
 * or taking the orbital term about a different origin than the spin term
 * flips one contribution and makes real errors "cancel" into a fake pass.
 * Every function here uses the convention above.
 */
#ifndef AM_H
#define AM_H

#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

#define AM_TOL_ONE_PER_MILLE 1e-3   /* requirement: |dL| / L_ref < 0.001 */

typedef struct { double x, y, z; } am_vec3;

am_vec3 am_v(double x, double y, double z);
am_vec3 am_add(am_vec3 a, am_vec3 b);
am_vec3 am_sub(am_vec3 a, am_vec3 b);
am_vec3 am_scale(am_vec3 a, double s);
am_vec3 am_cross(am_vec3 a, am_vec3 b);
double  am_dot(am_vec3 a, am_vec3 b);
double  am_len2(am_vec3 a);
double  am_len(am_vec3 a);

/* m <= 0 (or I <= 0) marks an infinite-mass / non-rotating static body. */
typedef struct {
    am_vec3 x;   /* COM world position [m] */
    am_vec3 v;   /* COM linear velocity [m/s] */
    am_vec3 w;   /* angular velocity [rad/s], right-hand signed */
    double  m;   /* mass [kg]; <= 0 -> static */
    double  I;   /* scalar (isotropic) inertia [kg m^2]; <= 0 -> no spin */
} am_body;

typedef struct {
    am_vec3 p;   /* shared world contact point [m] */
    am_vec3 n;   /* unit normal, A -> B */
    double  e;   /* restitution [0..1] */
} am_contact;

typedef enum {
    AM_TIER_COARSE = 0, /* particle solve: lever arm ignored, 0 cross products */
    AM_TIER_F32    = 1, /* full rotational arm, arithmetic in float32 */
    AM_TIER_F64    = 2  /* full rotational arm, arithmetic in double */
} am_tier;

typedef struct {
    int    status;   /* 1 impulse applied, 0 already separating, <0 bad args */
    double j;        /* normal impulse magnitude [N s] */
    double k_eff;    /* effective inverse mass used in j = -(1+e) vn / k */
    double vn;       /* pre-solve contact normal relative speed */
    am_vec3 rA, rB;  /* true lever arms COM -> p */
    double eccA;     /* eccentricity |rA x n| = perpendicular dist to n-line */
    double eccB;
} am_solution;

/* Snapshot of the pre-collision angular momentum about the pre-collision
 * system COM. The SAME origin must be used for before and after. */
typedef struct {
    am_vec3 origin;
    am_vec3 L_before;   /* signed vector total [kg m^2/s] */
    double  total_before; /* |L_before| */
    double  gross;      /* cancellation-free reference scale L_ref */
} am_snapshot;

typedef struct {
    am_vec3 L_before, L_after, dL;
    double  dL_abs;       /* |L_after - L_before| -- the "change amount" */
    double  total_before; /* |L_before| */
    double  gross;        /* L_ref */
    double  rel_total;    /* literal |dL| / |L_before| (HUGE if near cancel) */
    double  rel_gross;    /* enforced |dL| / L_ref */
    int     near_cancel;  /* |L_before| < 0.1 * gross */
    double  tol;
    double  eccA, eccB, ecc_max; /* eccentricities at the collision [m] */
    int     pass;         /* rel_gross <= tol */
} am_report;

am_snapshot am_snapshot_make(const am_body* A, const am_body* B);

/* Solve one frictionless normal impulse and apply it in place.
 * errA/errB are injected lever-arm errors added to the true arm for the
 * rotational update (pass am_v(0,0,0) for an exact arm). They model the
 * arm-computation precision under test. Ignored by AM_TIER_COARSE. */
am_solution am_solve_apply(am_body* A, am_body* B, const am_contact* c,
                           am_tier tier, am_vec3 errA, am_vec3 errB);

/* Independent conservation check from two full body states. */
am_report am_check(const am_snapshot* snap,
                   const am_body* A_after, const am_body* B_after,
                   double eccA, double eccB, double tol);

/* Total / gross angular momentum about an arbitrary fixed origin. */
am_vec3 am_total_L(const am_body* A, const am_body* B, am_vec3 origin);
double  am_gross_scale(const am_body* A, const am_body* B, am_vec3 origin);

/* Build the canonical worked collision used by the self-test and the
 * sensitivity calibration. Returns the pre-state and fills the contact.
 * Geometry is hand-checkable: L_before = (0,0,3.0), J ~= 5.970,
 * eccA = 0.30 m, eccB = 0.15 m. */
void am_canonical(am_body* A, am_body* B, am_contact* c);

#ifdef __cplusplus
}
#endif

#endif /* AM_H */
