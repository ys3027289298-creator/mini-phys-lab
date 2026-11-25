/*
 * degenerate: sphere-vs-sphere contact normal at centre coincidence.
 * See degenerate.h for the fallback contract.
 */
#include "degenerate.h"
#include <math.h>
#include <float.h>

float deg_pos_eps_ratio = 1e-6f;
float deg_vel_eps       = 1e-6f;
V3    deg_fixed_axis    = {1.0f, 0.0f, 0.0f};

static V3 v3_sub(V3 a, V3 b) {
    return (V3){a.x - b.x, a.y - b.y, a.z - b.z};
}
static float v3_dot(V3 a, V3 b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}
static float v3_len(V3 a) { return sqrtf(v3_dot(a, a)); }

static V3 v3_normalize_or(V3 v, V3 fallback) {
    float l = v3_len(v);
    if (!isfinite(l) || l < deg_vel_eps) return fallback;
    return (V3){v.x / l, v.y / l, v.z / l};
}

int v3_is_bad(V3 v) {
    return !isfinite(v.x) || !isfinite(v.y) || !isfinite(v.z);
}

int v3_is_unit(V3 v, float tol) {
    if (v3_is_bad(v)) return 0;
    return fabsf(v3_len(v) - 1.0f) <= tol;
}

void history_init(PairHistory *h, PairFrame *storage, int capacity) {
    h->ring = storage;
    h->capacity = capacity;
    h->head = 0;
    h->count = 0;
}

void history_push(PairHistory *h, V3 rel_pos, V3 rel_vel) {
    if (h->capacity <= 0) return;
    h->ring[h->head].rel_pos = rel_pos;
    h->ring[h->head].rel_vel = rel_vel;
    h->head = (h->head + 1) % h->capacity;
    if (h->count < h->capacity) h->count++;
}

/* Newest frame first: index 0 is the most recently pushed. */
static const PairFrame *history_at(const PairHistory *h, int age) {
    if (age >= h->count) return NULL;
    int idx = (h->head - 1 - age + h->capacity) % h->capacity;
    return &h->ring[idx];
}

/*
 * Why these fallbacks in this order:
 *
 * Let d(t) = cB(t) - cA(t). The contact normal is n = d / |d|.
 * If d -> 0 while d' = rel_vel stays finite and nonzero, l'Hopital gives
 * n -> d' / |d'|: the limit exists and is the heading of the relative
 * velocity, even though the position quotient is 0/0. That limit depends on
 * the *path* by which the pair reached coincidence, not on the position, so
 * it cannot be recovered from geometry alone.
 *
 * If d' -> 0 as well (the pair arrives dead still, or is born overlapping at
 * rest), both quotients vanish and no direction is derivable from the current
 * state at all. The only remaining physical evidence is where the pair was in
 * earlier frames: walk the history back until a frame carries a usable
 * position or velocity. When the whole window is empty, a fixed world axis is
 * the last resort -- deterministic and unit length, but geometrically
 * arbitrary, so the response pushes the spheres along the axis rather than
 * bouncing them apart.
 */
Contact sphere_contact(const Sphere *a, const Sphere *b,
                       V3 rel_vel, const PairHistory *h) {
    Contact out;
    out.normal = V3_ZERO;
    out.depth = 0.0f;
    out.hit = 0;
    out.source = NORMAL_FROM_POSITION;
    out.history_age = 0;

    V3 d = v3_sub(b->c, a->c);
    float dist = v3_len(d);
    float rsum = a->r + b->r;
    float pos_eps = rsum * deg_pos_eps_ratio;
    out.depth = rsum - dist;
    if (dist > rsum) return out;           /* separated: no contact */
    out.hit = 1;

    if (dist > pos_eps) {
        /* Ordinary case: geometry defines the normal, pointing B -> A. */
        out.normal = (V3){-d.x / dist, -d.y / dist, -d.z / dist};
        out.source = NORMAL_FROM_POSITION;
        return out;
    }

    /*
     * Coincident centres: the 0/0 branch. d -> 0, and with rel_vel = d'
     * finite and nonzero l'Hopital gives d/|d| -> d'/|d'|; the contact
     * normal is its opposite (B -> A convention).
     */
    if (v3_len(rel_vel) >= deg_vel_eps) {
        V3 h = v3_normalize_or(rel_vel, V3_ZERO);
        out.normal = (V3){-h.x, -h.y, -h.z};
        if (v3_is_unit(out.normal, 1e-6f)) {
            out.source = NORMAL_FROM_VELOCITY;
            return out;
        }
    }

    if (h) {
        for (int age = 0; age < h->count; age++) {
            const PairFrame *f = history_at(h, age);
            float plen = v3_len(f->rel_pos);
            if (plen > pos_eps) {
                out.normal = (V3){-f->rel_pos.x / plen,
                                  -f->rel_pos.y / plen,
                                  -f->rel_pos.z / plen};
                out.source = NORMAL_FROM_HISTORY;
                out.history_age = age + 1;
                return out;
            }
            float vlen = v3_len(f->rel_vel);
            if (vlen >= deg_vel_eps) {
                out.normal = (V3){-f->rel_vel.x / vlen,
                                  -f->rel_vel.y / vlen,
                                  -f->rel_vel.z / vlen};
                out.source = NORMAL_FROM_HISTORY;
                out.history_age = age + 1;
                return out;
            }
        }
    }

    /* No physical evidence anywhere in the window. */
    out.normal = v3_normalize_or(deg_fixed_axis, (V3){1, 0, 0});
    out.source = NORMAL_FROM_FIXED_AXIS;
    return out;
}

Contact sphere_contact_naive(const Sphere *a, const Sphere *b) {
    Contact out;
    V3 d = v3_sub(b->c, a->c);
    float dist = sqrtf(d.x*d.x + d.y*d.y + d.z*d.z);
    out.depth = a->r + b->r - dist;
    out.hit = dist <= a->r + b->r;
    /* Exact coincidence: 0/0 -> NaN normal, NaN impulse downstream. */
    out.normal = (V3){d.x / dist, d.y / dist, d.z / dist};
    out.source = NORMAL_FROM_POSITION;
    out.history_age = 0;
    return out;
}

void apply_impulse(V3 *va, V3 *vb, V3 n, float restitution) {
    /*
     * Equal masses, n points B -> A. vn = (vB - vA).n < 0 means B moves
     * against the push-out direction, i.e. the pair is approaching. Each
     * body takes half of the full relative-velocity reversal (1 + e)*|vn|.
     */
    V3 rv = v3_sub(*vb, *va);
    float vn = v3_dot(rv, n);
    if (vn >= 0.0f) return;              /* already separating */
    float j = (1.0f + restitution) * vn * 0.5f;   /* j < 0 */
    va->x += j * n.x; va->y += j * n.y; va->z += j * n.z;
    vb->x -= j * n.x; vb->y -= j * n.y; vb->z -= j * n.z;
}
