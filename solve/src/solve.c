/* solve: normal + tangential impulse resolution for a sphere hitting a plane.

   Body type is a required argument (see solve.h). Kinematic bodies have
   infinite mass: the solver must not modify their velocity. */
#include <stdio.h>
#include <math.h>

#include "solve.h"

/* Whether an impulse may change this body's velocity.
   Exhaustive over BodyKind: adding a kind without a policy is a -Wswitch error. */
static int takes_impulse(BodyType t) {
    switch (t.kind) {
    case BK_DYNAMIC:   return 1;
    case BK_KINEMATIC: return 0;
    }
    return 0;   /* unreachable; keeps -Wreturn-type quiet */
}

int resolve(BodyType t, V3* v, V3 n, const Material* m, float* out_normal_loss) {
    if (!takes_impulse(t)) {          /* kinematic: push others, never get pushed */
        if (out_normal_loss) *out_normal_loss = 0.0f;
        return 0;
    }
    if (m->restitution < 0.0f || m->restitution > 1.0f) return -1;
    if (m->friction < 0.0f || m->friction > 1.0f) return -1;

    float vn = v->x*n.x + v->y*n.y + v->z*n.z;
    if (vn > 0.0f) {          /* separating already */
        if (out_normal_loss) *out_normal_loss = 0.0f;
        return 0;
    }
    /* tangential part */
    V3 tng = { v->x - vn*n.x, v->y - vn*n.y, v->z - vn*n.z };
    float tl = sqrtf(tng.x*tng.x + tng.y*tng.y + tng.z*tng.z);

    /* normal impulse: reverse and scale */
    float new_vn = -vn * m->restitution;
    /* tangential impulse: reduce tangential speed by friction * |vn| */
    float reduce = m->friction * fabsf(vn);
    if (reduce > tl) reduce = tl;    /* friction cannot reverse tangential motion */
    float scale = (tl > 1e-9f) ? (tl - reduce) / tl : 0.0f;

    v->x = tng.x * scale + n.x * new_vn;
    v->y = tng.y * scale + n.y * new_vn;
    v->z = tng.z * scale + n.z * new_vn;

    if (out_normal_loss) *out_normal_loss = fabsf(vn) - fabsf(new_vn);
    return 1;
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
    check("resolve_ok", resolve(BODY_DYNAMIC, &v, n, &m, &loss) == 1);
    check("bounce_up", v.y > 0.0f);
    check("restitution_applied", fabsf(v.y - 5.0f) < 1e-3f);
    check("loss_recorded", fabsf(loss - 5.0f) < 1e-3f);

    /* separating velocity must be left alone */
    V3 up = { 0, 3, 0 };
    check("separating_untouched", resolve(BODY_DYNAMIC, &up, n, &m, &loss) == 0 && fabsf(up.y - 3.0f) < 1e-6f);

    /* friction reduces tangential, never reverses it */
    V3 slide = { 4, -1, 0 };
    Material grippy = { 0.0f, 1.0f };
    resolve(BODY_DYNAMIC, &slide, n, &grippy, &loss);
    check("tangential_not_reversed", slide.x >= 0.0f);
    check("tangential_reduced", slide.x < 4.0f);

    /* zero restitution: pure slide, no bounce */
    Material dull = { 0.0f, 0.0f };
    V3 flat = { 2, -2, 0 };
    resolve(BODY_DYNAMIC, &flat, n, &dull, &loss);
    check("no_bounce", fabsf(flat.y) < 1e-6f);

    /* regression: the incident itself. A kinematic body driven into the plane
       must come out bit-identical and untouched. */
    V3 kine = { 3, -7, 1 };
    float kine_loss = 123.0f;
    check("kinematic_untouched", resolve(BODY_KINEMATIC, &kine, n, &m, &kine_loss) == 0);
    check("kinematic_velocity_unchanged",
          kine.x == 3.0f && kine.y == -7.0f && kine.z == 1.0f);
    check("kinematic_no_impulse_loss", kine_loss == 0.0f);

    Material bad = { 1.5f, 0.2f };
    V3 any = { 0, -1, 0 };
    check("bad_restitution", resolve(BODY_DYNAMIC, &any, n, &bad, &loss) == -1);
    Material bad2 = { 0.5f, -0.1f };
    check("bad_friction", resolve(BODY_DYNAMIC, &any, n, &bad2, &loss) == -1);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
