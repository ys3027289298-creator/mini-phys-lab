
/* collide_ss 实现：见 collide_ss.h 中的判据与兜底说明。 */
#include "collide_ss.h"
#include <math.h>

static int valid_unit(const V3* n) {
    if (!isfinite(n->x) || !isfinite(n->y) || !isfinite(n->z)) return 0;
    float l2 = n->x * n->x + n->y * n->y + n->z * n->z;
    return fabsf(l2 - 1.0f) < 1e-5f;
}

int balls_contact(const Ball* a, const Ball* b, float eps_t, Contact* o) {
    if (a->r <= 0.0f || b->r <= 0.0f) return -1;
#ifndef SWEEP_ALLOW_ZERO_EPS
    if (eps_t <= 0.0f) eps_t = EPS_T_DEFAULT;
#endif

    float dx = b->c.x - a->c.x;
    float dy = b->c.y - a->c.y;
    float dz = b->c.z - a->c.z;
    float d2 = dx * dx + dy * dy + dz * dz;
    float rsum = a->r + b->r;

    o->hit = 0;
    o->depth = 0.0f;
    if (d2 >= rsum * rsum) return 0;

    /* 用 d2 与 eps_t^2 比较，退化路径上不做可能为零的 sqrt/除法 */
    if (d2 >= eps_t * eps_t) {
        float d = sqrtf(d2);
        o->normal.x = dx / d;
        o->normal.y = dy / d;
        o->normal.z = dz / d;
        o->mode = MODE_POSITION;
    } else {
        /* 一级兜底：相对速度。约定 n 由 A 指向 B，取 A 相对 B 的来流方向 */
        float rvx = a->v.x - b->v.x;
        float rvy = a->v.y - b->v.y;
        float rvz = a->v.z - b->v.z;
        float rv2 = rvx * rvx + rvy * rvy + rvz * rvz;
        if (rv2 > 0.0f && isfinite(rv2)) {
            float rvl = sqrtf(rv2);
            o->normal.x = rvx / rvl;
            o->normal.y = rvy / rvl;
            o->normal.z = rvz / rvl;
            o->mode = MODE_VELOCITY;
        } else {
            /* 二级兜底：两球同速（相对速度也退化），固定轴 */
            o->normal.x = 1.0f;
            o->normal.y = 0.0f;
            o->normal.z = 0.0f;
            o->mode = MODE_FIXED_AXIS;
        }
    }

    /* 不变量：任何接触都必须携带非零单位法线；违反则拒绝输出 */
    if (!valid_unit(&o->normal)) {
        o->hit = 0;
        o->depth = 0.0f;
        return -2;
    }

    o->depth = rsum - sqrtf(d2);   /* d2 可接受为 0：sqrt(0)=0，深度=rsum */
    o->hit = 1;
    return 1;
}
