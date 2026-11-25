
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

/* ---- angular-aware rolling contact (sphere vs fixed plane) ----
 * The point-mass resolve() above takes its tangent direction from v only, so
 * it cannot represent rolling slip u = v + omega x r and never updates omega.
 * resolve_rolling() is the additive path used by the rolling scenarios:
 * correct contact-point slip direction, rolling effective mass, angular
 * impulse r x J, and optional torsional friction about the normal. */
typedef struct {
    V3 v, w;        /* linear and angular velocity */
    float R;        /* sphere radius, > 0 */
    float inv_m;    /* 1/m */
    float inv_I;    /* 1/I (0 = point mass, no rotational response) */
} RollBody;

static int resolve_rolling(RollBody* b, V3 n, float friction, float mu_t) {
    if (friction < 0.0f || friction > 1.0f || mu_t < 0.0f) return -1;
    if (b->R <= 0.0f || b->inv_m <= 0.0f || b->inv_I < 0.0f) return -1;

    V3 r = { -b->R*n.x, -b->R*n.y, -b->R*n.z };       /* centre -> contact */
    V3 wxr = { b->w.y*r.z - b->w.z*r.y,
               b->w.z*r.x - b->w.x*r.z,
               b->w.x*r.y - b->w.y*r.x };
    V3 u = { b->v.x + wxr.x, b->v.y + wxr.y, b->v.z + wxr.z };
    float un = u.x*n.x + u.y*n.y + u.z*n.z;
    if (un > 0.0f) return 0;                           /* separating */

    /* normal impulse (e=0): k_n = inv_m since r x n = 0 */
    float jn = -un / b->inv_m;
    b->v.x += n.x * jn * b->inv_m;
    b->v.y += n.y * jn * b->inv_m;
    b->v.z += n.z * jn * b->inv_m;

    /* recompute contact slip after the normal impulse */
    wxr = (V3){ b->w.y*r.z - b->w.z*r.y,
                b->w.z*r.x - b->w.x*r.z,
                b->w.x*r.y - b->w.y*r.x };
    u = (V3){ b->v.x + wxr.x, b->v.y + wxr.y, b->v.z + wxr.z };
    float udotn = u.x*n.x + u.y*n.y + u.z*n.z;
    V3 ut = { u.x - udotn*n.x, u.y - udotn*n.y, u.z - udotn*n.z };
    float slip = sqrtf(ut.x*ut.x + ut.y*ut.y + ut.z*ut.z);

    if (slip > 1e-9f && jn > 0.0f) {
        /* tangent effective inverse mass k_t = inv_m + |r x t|^2 inv_I,
         * |r x t| = R for a tangent unit direction and isotropic inertia */
        float kt = b->inv_m + b->R*b->R*b->inv_I;
        float jmag = slip / kt;                        /* sticks the slip */
        float cap = friction * jn;                     /* Coulomb cone */
        if (jmag > cap) jmag = cap;
        float s = -jmag / slip;                        /* OPPOSE slip */
        V3 J = { ut.x*s, ut.y*s, ut.z*s };
        b->v.x += J.x * b->inv_m;
        b->v.y += J.y * b->inv_m;
        b->v.z += J.z * b->inv_m;
        V3 rxJ = { r.y*J.z - r.z*J.y,
                   r.z*J.x - r.x*J.z,
                   r.x*J.y - r.y*J.x };
        b->w.x += rxJ.x * b->inv_I;
        b->w.y += rxJ.y * b->inv_I;
        b->w.z += rxJ.z * b->inv_I;
    }

    /* torsional friction: damp spin component along n, capped by mu_t R jn */
    if (mu_t > 0.0f && b->inv_I > 0.0f) {
        float wn = b->w.x*n.x + b->w.y*n.y + b->w.z*n.z;
        float I = 1.0f / b->inv_I;
        float H = fminf(fabsf(I*wn), mu_t*b->R*jn);
        float dwn = copysignf(H*b->inv_I, wn);
        b->w.x -= n.x*dwn;
        b->w.y -= n.y*dwn;
        b->w.z -= n.z*dwn;
    }
    return 1;
}

/* magnitude of tangential contact slip for the tests below */
static float rolling_slip(const RollBody* b, V3 n) {
    V3 r = { -b->R*n.x, -b->R*n.y, -b->R*n.z };
    V3 wxr = { b->w.y*r.z - b->w.z*r.y,
               b->w.z*r.x - b->w.x*r.z,
               b->w.x*r.y - b->w.y*r.x };
    V3 u = { b->v.x + wxr.x, b->v.y + wxr.y, b->v.z + wxr.z };
    float un = u.x*n.x + u.y*n.y + u.z*n.z;
    V3 ut = { u.x-un*n.x, u.y-un*n.y, u.z-un*n.z };
    return sqrtf(ut.x*ut.x + ut.y*ut.y + ut.z*ut.z);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
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

    /* ---- rolling resolver checks (solid sphere R=.5, m=1, I=.1) ---- */
    RollBody rb; memset(&rb, 0, sizeof rb);
    rb.R = 0.5f; rb.inv_m = 1.0f; rb.inv_I = 10.0f;
    rb.v = (V3){ 2.0f, 0.0f, 0.0f };
    rb.w = (V3){ 0.0f, 0.0f, -4.0f };       /* pure roll: v + w x r = 0 */
    check("pure_roll_no_slip", rolling_slip(&rb, n) < 1e-5f);
    resolve_rolling(&rb, n, 0.8f, 0.0f);
    check("pure_roll_untouched", fabsf(rb.v.x - 2.0f) < 1e-4f);

    /* wobble creates lateral slip; impulse must shrink it and update omega */
    RollBody wob; memset(&wob, 0, sizeof wob);
    wob.R = 0.5f; wob.inv_m = 1.0f; wob.inv_I = 10.0f;
    wob.v = (V3){ 2.0f, -1.0f, 0.0f };
    wob.w = (V3){ 0.5f, 0.0f, -4.0f };
    float slip0 = rolling_slip(&wob, n);
    float wx0 = wob.w.x;
    resolve_rolling(&wob, n, 0.8f, 0.0f);
    check("wobble_slip_reduced", rolling_slip(&wob, n) < slip0);
    check("wobble_omega_updated", fabsf(wob.w.x) < fabsf(wx0));
    check("wobble_lateral_impulse_dir", wob.v.z > 0.0f);

    /* Coulomb cap: low friction trims slip but never reverses it */
    RollBody cap; memset(&cap, 0, sizeof cap);
    cap.R = 0.5f; cap.inv_m = 1.0f; cap.inv_I = 10.0f;
    cap.v = (V3){ 0.0f, -2.0f, 5.0f };
    float sz_before = 5.0f;
    resolve_rolling(&cap, n, 0.05f, 0.0f);
    check("friction_cone_no_reverse", cap.v.z >= 0.0f && cap.v.z < sz_before);

    /* torsional friction kills normal-axis spin */
    RollBody spin; memset(&spin, 0, sizeof spin);
    spin.R = 0.5f; spin.inv_m = 1.0f; spin.inv_I = 10.0f;
    spin.v = (V3){ 0.0f, -1.0f, 0.0f };
    spin.w = (V3){ 0.0f, 6.0f, 0.0f };
    resolve_rolling(&spin, n, 0.8f, 0.5f);
    check("torsion_spin_damped", fabsf(spin.w.y) < 6.0f);

    check("rolling_bad_args", resolve_rolling(&rb, n, 1.5f, 0.0f) == -1);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
