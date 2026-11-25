
/* solve: normal + tangential impulse resolution for a sphere hitting a plane.
 *
 * Friction model: Coulomb friction cone.  The normal impulse is
 * jn = (1+e)|vn| (unit mass), and the cone bounds the tangential impulse by
 *
 *     |jt| <= mu * jn,        tan(alpha) = mu
 *
 * where alpha is the cone half-angle (full opening angle = 2*alpha).
 * Each tangent axis is clipped independently (the pyramid form of the cone):
 * the same absolute amount delta = mu*(1+e)*|vn| is removed from every
 * tangential component, clamped so friction never reverses a component.
 * Therefore the *fraction* removed from component i is
 *
 *     r_i = min(1, mu*(1+e)*|vn| / |vt_i|)
 *
 * which is direction-dependent: small tangential components stick first.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

typedef struct { float x, y, z; } V3;

typedef struct {
    float restitution;   /* 0..1 */
    float friction;      /* mu, 0..1 */
} Material;

static V3 v3cross(V3 a, V3 b) {
    V3 r = { a.y*b.z - a.z*b.y,
             a.z*b.x - a.x*b.z,
             a.x*b.y - a.y*b.x };
    return r;
}

static float v3dot(V3 a, V3 b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

static V3 v3scale(V3 a, float s) {
    V3 r = { a.x*s, a.y*s, a.z*s };
    return r;
}

static V3 v3norm(V3 a) {
    float l = sqrtf(v3dot(a, a));
    return l > 1e-12f ? v3scale(a, 1.0f / l) : (V3){ 0, 0, 0 };
}

/* orthonormal tangent basis for a unit normal n */
static void tangent_basis(V3 n, V3* t1, V3* t2) {
    V3 ref = fabsf(n.y) > 0.9f ? (V3){ 1, 0, 0 } : (V3){ 0, 1, 0 };
    *t1 = v3norm(v3cross(ref, n));
    *t2 = v3cross(n, *t1);
}

/* cone geometry: mu -> half-angle and full opening angle, degrees.
 * fails when mu is negative or so large the angle degenerates to 90 deg. */
static int cone_angles(float mu, float* half_deg, float* opening_deg) {
    if (mu < 0.0f || !isfinite(mu)) return -3;
    float half = atanf(mu) * 180.0f / 3.14159265358979f;
    if (half >= 90.0f) return -3;       /* tan(90) = inf: unbounded tangential impulse */
    if (half_deg) *half_deg = half;
    if (opening_deg) *opening_deg = 2.0f * half;
    return 0;
}

/* inverse design: pick mu from a wanted tangential reduction ratio r at a
 * known incidence.  mu = r*|vt| / ((1+e)*|vn|), then alpha = atan(mu). */
static int design_friction(float target_r, float vt, float vn,
                           float restitution, float* mu_out,
                           float* half_deg, float* opening_deg) {
    if (target_r <= 0.0f || target_r > 1.0f) return -4;
    if (vt <= 0.0f || vn >= 0.0f) return -4;
    if (restitution < 0.0f || restitution > 1.0f) return -1;
    float mu = target_r * vt / ((1.0f + restitution) * (-vn));
    int rc = cone_angles(mu, half_deg, opening_deg);
    if (rc != 0) return rc;
    if (mu > 1.0f) return -3;           /* material coefficient out of range */
    if (mu_out) *mu_out = mu;
    return 0;
}

/* build a material from a chosen cone half-angle; rejects alpha >= 90 deg */
static int material_from_cone(float half_deg, float restitution, Material* out) {
    if (restitution < 0.0f || restitution > 1.0f) return -1;
    if (!(half_deg >= 0.0f && half_deg < 90.0f)) return -3;
    float rad = half_deg * 3.14159265358979f / 180.0f;
    float mu = tanf(rad);
    if (!isfinite(mu) || mu < 0.0f || mu > 1.0f) return -3;
    out->restitution = restitution;
    out->friction = mu;
    return 0;
}

/* resolve one contact against a unit normal pointing out of the surface.
 * returns 1 = impulse applied, 0 = already separating, <0 = bad material. */
static int resolve(V3* v, V3 n, const Material* m, float* out_normal_loss) {
    if (m->restitution < 0.0f || m->restitution > 1.0f) return -1;
    if (m->friction < 0.0f || m->friction > 1.0f) return -1;

    float vn = v3dot(*v, n);
    if (vn > 0.0f) {          /* separating already */
        if (out_normal_loss) *out_normal_loss = 0.0f;
        return 0;
    }

    /* normal impulse reverses and scales the normal component */
    float new_vn = -vn * m->restitution;

    /* friction cone: per-axis tangential impulse cap mu * jn */
    float delta = m->friction * (1.0f + m->restitution) * fabsf(vn);

    V3 t1, t2;
    tangent_basis(n, &t1, &t2);
    float a = v3dot(*v, t1);
    float b = v3dot(*v, t2);

    /* clip each tangential component independently, never reversing it */
    float a2 = fabsf(a) > delta ? a - copysignf(delta, a) : 0.0f;
    float b2 = fabsf(b) > delta ? b - copysignf(delta, b) : 0.0f;

    v->x = t1.x * a2 + t2.x * b2 + n.x * new_vn;
    v->y = t1.y * a2 + t2.y * b2 + n.y * new_vn;
    v->z = t1.z * a2 + t2.z * b2 + n.z * new_vn;

    if (out_normal_loss) *out_normal_loss = fabsf(vn) - fabsf(new_vn);
    return 1;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* theoretical per-component reduction ratio, including the stick cap */
static float theory_ratio(float vti, float vn, const Material* m) {
    float r = m->friction * (1.0f + m->restitution) * fabsf(vn) / fabsf(vti);
    return r > 1.0f ? 1.0f : r;
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

    /* ---- cone geometry: mu = 0.4 -> half angle atan(0.4) -------------- */
    Material cm = { 0.3f, 0.4f };
    float half = 0, opening = 0;
    check("cone_angles_ok", cone_angles(cm.friction, &half, &opening) == 0);
    printf("cone half=%.4f deg opening=%.4f deg\n", half, opening);
    check("half_angle_atan_mu", fabsf(half - atanf(0.4f) * 180.0f / 3.14159265f) < 1e-3f);
    check("opening_is_double", fabsf(opening - 2.0f * half) < 1e-3f);

    /* inverse design round trip: pick r at 30 deg incidence, recover alpha */
    float theta = 30.0f * 3.14159265f / 180.0f;
    float speed = 10.0f;
    float vn30 = -speed * sinf(theta);
    float vt30 = speed * cosf(theta);
    float target_r = theory_ratio(vt30, vn30, &cm);
    float mu_d = 0, half_d = 0, open_d = 0;
    check("design_invert_ok",
          design_friction(target_r, vt30, vn30, cm.restitution,
                          &mu_d, &half_d, &open_d) == 0);
    printf("design target_r=%.4f -> mu=%.4f half=%.4f deg\n", target_r, mu_d, half_d);
    check("design_recovers_mu", fabsf(mu_d - cm.friction) < 1e-3f);
    check("design_recovers_angle", fabsf(half_d - half) < 1e-2f);

    /* material_from_cone is the exact inverse of cone_angles */
    Material rebuilt;
    check("from_cone_ok", material_from_cone(half, cm.restitution, &rebuilt) == 0);
    check("from_cone_mu", fabsf(rebuilt.friction - cm.friction) < 1e-4f);

    /* ---- three incidence angles, two orthogonal tangent components ----- */
    V3 t1, t2;
    tangent_basis(n, &t1, &t2);
    float phi = 35.0f * 3.14159265f / 180.0f;   /* split of tangential speed */
    float angles[3] = { 15.0f, 30.0f, 45.0f };
    int angle_ok = 1;
    for (int k = 0; k < 3; k++) {
        float th = angles[k] * 3.14159265f / 180.0f;
        float vnv = -speed * sinf(th);
        float vtv = speed * cosf(th);
        float a0 = vtv * cosf(phi);
        float b0 = vtv * sinf(phi);
        V3 shot = {
            n.x * vnv + t1.x * a0 + t2.x * b0,
            n.y * vnv + t1.y * a0 + t2.y * b0,
            n.z * vnv + t1.z * a0 + t2.z * b0
        };
        resolve(&shot, n, &cm, &loss);
        float a1 = v3dot(shot, t1);
        float b1 = v3dot(shot, t2);
        float r1 = (a0 - a1) / a0;             /* measured reduction ratios */
        float r2 = (b0 - b1) / b0;
        float e1 = theory_ratio(a0, vnv, &cm);
        float e2 = theory_ratio(b0, vnv, &cm);
        /* 3% relative tolerance; absolute 0.03 where theory hits the cap */
        float tol1 = e1 >= 1.0f ? 0.03f : 0.03f * e1;
        float tol2 = e2 >= 1.0f ? 0.03f : 0.03f * e2;
        int ok1 = fabsf(r1 - e1) <= tol1;
        int ok2 = fabsf(r2 - e2) <= tol2;
        printf("theta=%2.0f r1_meas=%.4f r1_theory=%.4f r2_meas=%.4f r2_theory=%.4f\n",
               angles[k], r1, e1, r2, e2);
        if (!ok1 || !ok2) {
            angle_ok = 0;
            printf("  MISMATCH actual r1=%.4f (want %.4f) r2=%.4f (want %.4f)\n",
                   r1, e1, r2, e2);
        }
    }
    check("three_angles_within_3pct", angle_ok);

    /* steeper shot: the small component must stick (cap at r = 1) while the
     * large one keeps sliding, proving the response is directional */
    float th = 55.0f * 3.14159265f / 180.0f;
    float vnv = -speed * sinf(th), vtv = speed * cosf(th);
    float a0 = vtv * cosf(phi), b0 = vtv * sinf(phi);
    V3 steep = {
        n.x * vnv + t1.x * a0 + t2.x * b0,
        n.y * vnv + t1.y * a0 + t2.y * b0,
        n.z * vnv + t1.z * a0 + t2.z * b0
    };
    resolve(&steep, n, &cm, &loss);
    float r1s = (a0 - v3dot(steep, t1)) / a0;
    float r2s = (b0 - v3dot(steep, t2)) / b0;
    printf("theta=55 r1=%.4f r2=%.4f (small component sticks first)\n", r1s, r2s);
    check("directional_stick", r2s >= 1.0f - 1e-5f && r1s < 0.95f);

    /* ---- infeasible parameter: 90 deg cone means mu is infinite -------- */
    Material badcone;
    check("cone90_rejected", material_from_cone(90.0f, 0.3f, &badcone) == -3);
    check("cone_over90_rejected", material_from_cone(95.0f, 0.3f, &badcone) == -3);
    check("cone_negative_rejected", material_from_cone(-1.0f, 0.3f, &badcone) == -3);
    check("inf_mu_rejected", cone_angles(INFINITY, &half, &opening) == -3);

    /* physical degeneracy at the limit: e = 0, mu -> huge kills every
     * component, so the ball leaves the wall with ~zero velocity */
    Material sticky = { 0.0f, 1.0f - 1e-6f };
    cone_angles(sticky.friction, &half, &opening);
    V3 cling = { 3.0f, -4.0f, 2.0f };
    resolve(&cling, n, &sticky, &loss);
    float residual = sqrtf(v3dot(cling, cling));
    printf("limit half=%.4f deg residual_speed=%.6f (sticks to wall)\n", half, residual);
    check("limit_sticks", residual < 1e-3f);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
