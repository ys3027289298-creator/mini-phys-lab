
/* solve: normal + tangential impulse resolution for a sphere hitting a plane. */
#include <stdio.h>
#include <string.h>
#include <math.h>

typedef struct { float x, y, z; } V3;

typedef struct {
    float restitution;   /* 0..1 */
    float friction;      /* 0..1 */
} Material;

/* resolve one contact against a unit normal pointing out of the surface */
static int resolve(V3* v, V3 n, const Material* m, float* out_normal_loss) {
    if (m->restitution < 0.0f || m->restitution > 1.0f) return -1;
    if (m->friction < 0.0f || m->friction > 1.0f) return -1;

    float vn = v->x*n.x + v->y*n.y + v->z*n.z;
    if (vn > 0.0f) {          /* separating already */
        if (out_normal_loss) *out_normal_loss = 0.0f;
        return 0;
    }
    /* tangential part */
    V3 t = { v->x - vn*n.x, v->y - vn*n.y, v->z - vn*n.z };
    float tl = sqrtf(t.x*t.x + t.y*t.y + t.z*t.z);

    /* normal impulse: reverse and scale */
    float new_vn = -vn * m->restitution;
    /* tangential impulse: reduce tangential speed by friction * |vn| */
    float reduce = m->friction * fabsf(vn);
    if (reduce > tl) reduce = tl;    /* friction cannot reverse tangential motion */
    float scale = (tl > 1e-9f) ? (tl - reduce) / tl : 0.0f;

    v->x = t.x * scale + n.x * new_vn;
    v->y = t.y * scale + n.y * new_vn;
    v->z = t.z * scale + n.z * new_vn;

    if (out_normal_loss) *out_normal_loss = fabsf(vn) - fabsf(new_vn);
    return 1;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* ---- 3D sphere/angular-contact regression ---- */

typedef struct { double x, y, z; } D3;

static D3 d3(double x, double y, double z) { D3 r = {x,y,z}; return r; }
static D3 dadd(D3 a, D3 b) { return d3(a.x+b.x, a.y+b.y, a.z+b.z); }
static D3 dsub(D3 a, D3 b) { return d3(a.x-b.x, a.y-b.y, a.z-b.z); }
static D3 dscale(D3 a, double s) { return d3(a.x*s, a.y*s, a.z*s); }
static double ddot(D3 a, D3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static D3 dcross(D3 a, D3 b) {
    return d3(a.y*b.z - a.z*b.y,
              a.z*b.x - a.x*b.z,
              a.x*b.y - a.y*b.x);
}
static double dlen(D3 a) { return sqrt(ddot(a,a)); }

typedef struct {
    int id;
    double radius;
    double inv_mass;
    double inv_inertia; /* isotropic sphere inertia inverse */
    D3 p, v, w;
    int awake;
} SpinBall;

typedef struct {
    D3 point;
    D3 normal;
    D3 tangent;
    int has_tangent;
} SpinContact;

enum {
    SPIN_APPLIED = 0,
    SPIN_BAD_MATERIAL = -1,
    SPIN_BAD_NORMAL = -2,
    SPIN_BAD_EFFECTIVE_MASS = -3,
    SPIN_SEPARATING = 1,
    SPIN_SLEEPING = 2,
    SPIN_ZERO_ARM_GUARD = 3
};

typedef struct {
    int reason;
    int normal_flipped;
    double vn_before;
    double vt_before;
    double jn, jt;
    double arm;
    double arm_relative_error;
    D3 r;
    D3 w_before, w_after;
} SpinReport;

typedef struct {
    double linear_speed_eps;
    double angular_speed_eps;
    double delay_s;
    double accumulated_s;
} SleepPolicy;

static void wake_ball(SpinBall* b, SleepPolicy* policy) {
    b->awake = 1;
    if (policy) policy->accumulated_s = 0.0;
}

/* Fixed policy: zero spin alone is not enough; translation must also be low. */
static int sleep_policy_tick(SpinBall* b, double dt, SleepPolicy* policy) {
    if (dlen(b->v) <= policy->linear_speed_eps &&
        dlen(b->w) <= policy->angular_speed_eps) {
        policy->accumulated_s += dt;
    } else {
        policy->accumulated_s = 0.0;
    }
    if (policy->accumulated_s >= policy->delay_s) b->awake = 0;
    return b->awake;
}

/* Some legacy paths guard a degenerate arm by zeroing the whole impulse. */
static int legacy_zero_arm_guard(D3 r, D3 tangent, double* j_out) {
    if (dlen(dcross(r, tangent)) < 1e-9) {
        *j_out = 0.0;
        return SPIN_ZERO_ARM_GUARD;
    }
    *j_out = 1.0; /* arm is usable; exact magnitude is irrelevant to the branch */
    return SPIN_APPLIED;
}

static const char* exclusion_step(const SpinReport* report) {
    switch (report->reason) {
        case SPIN_SLEEPING:
            return "sleep_gate";
        case SPIN_SEPARATING:
            return "vn_dot_sign_gate";
        case SPIN_ZERO_ARM_GUARD:
            return "lever_arm_guard";
        case SPIN_BAD_NORMAL:
            return "normal_unit_check";
        case SPIN_BAD_EFFECTIVE_MASS:
            return "effective_mass_denominator";
        case SPIN_APPLIED:
            if (report->arm < 1e-9) return "r_cross_torque_arm";
            if (dlen(dsub(report->w_after, report->w_before)) < 1e-12)
                return "angular_impulse_integration";
            return "none";
        default:
            return "parameter_validation";
    }
}

/*
 * Correct impulse-based sphere contact.  The contact vector r is from centre
 * to contact point.  A tangential impulse changes spin by invI * (r x J).
 */
static int solve_sphere_spin_contact(SpinBall* b, SpinContact* contact,
                                     double restitution, double friction,
                                     SpinReport* report) {
    memset(report, 0, sizeof(*report));
    report->reason = SPIN_APPLIED;
    report->w_before = b->w;

    if (restitution < 0.0 || restitution > 1.0 ||
        friction < 0.0 || friction > 1.0) {
        report->reason = SPIN_BAD_MATERIAL;
        return -1;
    }
    if (!b->awake) {
        report->reason = SPIN_SLEEPING;
        return 0;
    }

    double normal_len = dlen(contact->normal);
    if (!(normal_len > 1e-12)) {
        report->reason = SPIN_BAD_NORMAL;
        return -1;
    }
    contact->normal = dscale(contact->normal, 1.0 / normal_len);

    D3 centre_offset = dsub(b->p, contact->point);
    if (ddot(centre_offset, contact->normal) < 0.0) {
        contact->normal = dscale(contact->normal, -1.0);
        report->normal_flipped = 1;
    }
    D3 n = contact->normal;

    D3 r = dsub(contact->point, b->p);
    report->r = r;

    D3 contact_velocity = dadd(b->v, dcross(b->w, r));
    double vn = ddot(contact_velocity, n);
    report->vn_before = vn;
    if (vn > 0.0) {
        report->reason = SPIN_SEPARATING;
        return 0;
    }

    D3 tangent_velocity = dsub(contact_velocity, dscale(n, vn));
    double vt = dlen(tangent_velocity);
    report->vt_before = vt;

    D3 t = d3(0,0,0);
    if (vt > 1e-9) {
        t = dscale(tangent_velocity, 1.0 / vt);
    } else if (contact->has_tangent) {
        double tl = dlen(contact->tangent);
        if (tl > 1e-9) {
            t = dscale(contact->tangent, 1.0 / tl);
            t = dsub(t, dscale(n, ddot(t, n)));
            tl = dlen(t);
            if (tl > 1e-9) t = dscale(t, 1.0 / tl);
            else t = d3(0,0,0);
        }
    }

    D3 rn = dcross(r, n);
    double normal_effective_mass =
        b->inv_mass + b->inv_inertia * ddot(rn, rn);
    if (!(normal_effective_mass > 0.0)) {
        report->reason = SPIN_BAD_EFFECTIVE_MASS;
        return -1;
    }

    double jn = -(1.0 + restitution) * vn / normal_effective_mass;
    if (jn < 0.0 && jn > -1e-12) jn = 0.0;
    report->jn = jn;

    double jt = 0.0;
    if (dlen(t) > 0.5) {
        D3 rt = dcross(r, t);
        report->arm = dlen(rt);
        double tangent_effective_mass =
            b->inv_mass + b->inv_inertia * ddot(rt, rt);
        double required = -ddot(contact_velocity, t) / tangent_effective_mass;
        double limit = friction * jn;
        if (fabs(required) <= limit) {
            jt = required;
        } else {
            jt = copysign(limit, required);
        }
    }
    report->jt = jt;

    D3 impulse = dadd(dscale(n, jn), dscale(t, jt));
    b->v = dadd(b->v, dscale(impulse, b->inv_mass));

    D3 angular_impulse = dcross(r, impulse);
    b->w = dadd(b->w, dscale(angular_impulse, b->inv_inertia));
    report->w_after = b->w;

    D3 delta_w = dsub(report->w_after, report->w_before);
    D3 tangential_impulse = dscale(t, jt);
    double analytic_delta_w =
        dlen(dcross(r, tangential_impulse)) * b->inv_inertia;
    double actual_delta_w = dlen(delta_w);
    if (analytic_delta_w > 1e-12) {
        report->arm_relative_error =
            fabs(actual_delta_w - analytic_delta_w) / analytic_delta_w;
    } else {
        report->arm_relative_error = actual_delta_w <= 1e-12 ? 0.0 : 1.0;
    }

    report->reason = SPIN_APPLIED;
    return 1;
}

static SpinBall angular_test_ball(int id, int awake) {
    SpinBall b;
    memset(&b, 0, sizeof(b));
    b.id = id;
    b.radius = 0.5;
    b.inv_mass = 1.0;
    /* Solid sphere: I = 2/5 m R^2 = 0.1, hence invI = 10 for m=1,R=.5. */
    b.inv_inertia = 10.0;
    b.p = d3(0.0, 0.5, 0.0);
    b.v = d3(3.0, -1.0, 0.0);
    b.awake = awake;
    return b;
}

static void check_angular_fix(int ball_id, const char* cause,
                              const SpinBall* b, const SpinReport* report,
                              double expected_wz) {
    double rolling_ratio = fabs(b->w.z) * b->radius / 3.0;
    double ratio_error = fabs(b->w.z - expected_wz) / fabs(expected_wz);
    int fixed = report->reason == SPIN_APPLIED &&
                fabs(b->w.z) > 1e-6 &&
                ratio_error <= 0.05 &&
                report->arm_relative_error <= 0.05;

    check(ball_id == 1 ? "angular_fix_ball_1_zero_arm" :
          ball_id == 2 ? "angular_fix_ball_2_reversed_normal" :
                         "angular_fix_ball_3_sleep", fixed);
    printf("angular_ball=%d cause=%s wz=%+.6f analytic_wz=%+.6f "
           "rolling_ratio=%.5f analytic_ratio=0.71429 ratio_error=%.4f "
           "arm_relerr=%.2e status=%s\n",
           ball_id, cause, b->w.z, expected_wz, rolling_ratio,
           ratio_error, report->arm_relative_error,
           fixed ? "fixed" : "UNRESOLVED");
    if (!fixed) {
        printf("unresolved_ball=%d excluded_through=%s\n",
               ball_id, exclusion_step(report));
    }
}

static void run_angular_regression(void) {
    const double radius = 0.5;
    const double speed_x = 3.0;
    /* Sticking friction for a solid sphere gives omega_z = -5 U / (7 R). */
    const double expected_wz = -5.0 * speed_x / (7.0 * radius);

    /* Ball 1: contact point was stored at the centre, zeroing r x t. */
    SpinBall b1 = angular_test_ball(1, 1);
    SpinContact c1;
    memset(&c1, 0, sizeof(c1));
    c1.normal = d3(0.0, 1.0, 0.0);
    c1.point = b1.p; /* wrong: lever arm is zero */
    double legacy_j1 = -1.0;
    int legacy1 = legacy_zero_arm_guard(dsub(c1.point, b1.p),
                                        d3(1,0,0), &legacy_j1);
    check("zero_arm_legacy_path_sets_j_zero",
          legacy1 == SPIN_ZERO_ARM_GUARD && legacy_j1 == 0.0);

    c1.point = dadd(b1.p, dscale(c1.normal, -radius)); /* repaired surface point */
    SpinReport r1;
    solve_sphere_spin_contact(&b1, &c1, 0.0, 1.0, &r1);
    check("zero_arm_contact_repaired_to_surface", fabs(r1.arm - radius) < 1e-9);
    check_angular_fix(1, "zero_lever_arm", &b1, &r1, expected_wz);

    /* Ball 2: solver-side normal was flipped, making vn positive. */
    V3 legacy_v = {3.0f, -1.0f, 0.0f};
    V3 reversed_n = {0.0f, -1.0f, 0.0f};
    Material legacy_material = {0.0f, 1.0f};
    check("reversed_normal_legacy_separating_branch",
          resolve(&legacy_v, reversed_n, &legacy_material, NULL) == 0);

    SpinBall b2 = angular_test_ball(2, 1);
    SpinContact c2;
    memset(&c2, 0, sizeof(c2));
    c2.normal = d3(0.0, -1.0, 0.0); /* wrong input */
    c2.point = d3(0.0, 0.0, 0.0);
    SpinReport r2;
    solve_sphere_spin_contact(&b2, &c2, 0.0, 1.0, &r2);
    check("reversed_normal_is_reoriented", r2.normal_flipped == 1);
    check_angular_fix(2, "reversed_normal_positive_vn", &b2, &r2, expected_wz);

    /* Ball 3: over-eager sleep test used w == 0 alone and skipped solving. */
    SpinBall b3 = angular_test_ball(3, 0);
    SpinContact c3;
    memset(&c3, 0, sizeof(c3));
    c3.normal = d3(0.0, 1.0, 0.0);
    c3.point = d3(0.0, 0.0, 0.0);
    SpinReport r3;
    int first3 = solve_sphere_spin_contact(&b3, &c3, 0.0, 1.0, &r3);
    check("sleeping_ball_skipped",
          first3 == 0 && r3.reason == SPIN_SLEEPING &&
          r3.jn == 0.0 && r3.jt == 0.0);

    SleepPolicy policy = {0.05, 0.05, 0.5, 0.5};
    wake_ball(&b3, &policy); /* contact events wake/restart the body */
    int awake = sleep_policy_tick(&b3, 0.5, &policy);
    check("translating_zero_spin_ball_not_sleepy", awake == 1);

    memset(&c3, 0, sizeof(c3));
    c3.normal = d3(0.0, 1.0, 0.0);
    c3.point = d3(0.0, 0.0, 0.0);
    solve_sphere_spin_contact(&b3, &c3, 0.0, 1.0, &r3);
    check_angular_fix(3, "false_sleep_skip", &b3, &r3, expected_wz);
}

int main(void) {
    Material m = { 0.5f, 0.2f };
    V3 n = { 0, 1, 0 };

    V3 v = { 0, -10, 0 };
    float loss = 0;
    check("resolve_ok", resolve(&v, n, &m, &loss) == 1);
    check("bounce_up", v.y > 0.0f);
    check("restitution_applied", fabsf(v.y - 5.0f) < 1e-3f);
    check("loss_recorded", fabsf(loss - 5.0f) < 1e-3f);

    /* separating velocity must be left alone */
    V3 up = { 0, 3, 0 };
    check("separating_untouched", resolve(&up, n, &m, &loss) == 0 && fabsf(up.y - 3.0f) < 1e-6f);

    /* friction reduces tangential, never reverses it */
    V3 slide = { 4, -1, 0 };
    Material grippy = { 0.0f, 1.0f };
    resolve(&slide, n, &grippy, &loss);
    check("tangential_not_reversed", slide.x >= 0.0f);
    check("tangential_reduced", slide.x < 4.0f);

    /* zero restitution: pure slide, no bounce */
    Material dull = { 0.0f, 0.0f };
    V3 flat = { 2, -2, 0 };
    resolve(&flat, n, &dull, &loss);
    check("no_bounce", fabsf(flat.y) < 1e-6f);

    Material bad = { 1.5f, 0.2f };
    V3 any = { 0, -1, 0 };
    check("bad_restitution", resolve(&any, n, &bad, &loss) == -1);
    Material bad2 = { 0.5f, -0.1f };
    check("bad_friction", resolve(&any, n, &bad2, &loss) == -1);

    run_angular_regression();

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
