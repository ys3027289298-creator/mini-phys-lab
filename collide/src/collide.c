
/* collide: sphere vs AABB contact generation (normal + depth). */
#include <stdio.h>
#include <string.h>
#include <math.h>

typedef struct { float x, y, z; } V3;
typedef struct { V3 c; float r; } Sphere;
typedef struct { V3 min, max; } Aabb;

typedef struct { V3 normal; float depth; V3 point; int hit; } Contact;

typedef struct {
    float normal_length;
    float center_distance;
    float center_distance_from_depth;
    float depth_from_state;
    float depth_error;
    int normal_unit;
    int orientation_ok;
    int depth_ok;
    int ok;
} DepthAudit;

static float vdot(V3 a, V3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static V3 vsub(V3 a, V3 b) {
    V3 r = { a.x - b.x, a.y - b.y, a.z - b.z };
    return r;
}

static float vnorm(V3 a) {
    return sqrtf(vdot(a, a));
}

/*
 * Check a logged plane/sphere contact.  plane_point may be any point on the
 * plane.  For an outward unit normal n, h = dot(center - plane_point, n),
 * geometric depth = radius - h, and the log implies h = radius - depth.
 */
static DepthAudit audit_plane_depth(V3 center, V3 plane_point, float radius,
                                    Contact logged, float depth_tolerance) {
    DepthAudit a;
    memset(&a, 0, sizeof(a));

    float n_len = vnorm(logged.normal);
    a.normal_length = n_len;
    a.normal_unit = fabsf(n_len - 1.0f) <= 1e-4f;

    V3 n = logged.normal;
    if (n_len > 1e-12f) {
        n.x /= n_len; n.y /= n_len; n.z /= n_len;
    }

    float h = vdot(vsub(center, plane_point), n);
    a.center_distance = h;
    a.center_distance_from_depth = radius - logged.depth;
    a.depth_from_state = radius - h;
    a.depth_error = a.depth_from_state - logged.depth;

    /* Outward normal must point from the plane toward the sphere centre. */
    a.orientation_ok = h >= -depth_tolerance;
    a.depth_ok = fabsf(a.depth_error) <= depth_tolerance;
    a.ok = a.normal_unit && a.orientation_ok && a.depth_ok;
    return a;
}

static float clampf(float v, float lo, float hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

/* nearest point on box to sphere centre */
static V3 closest(const Aabb* b, V3 p) {
    V3 q;
    q.x = clampf(p.x, b->min.x, b->max.x);
    q.y = clampf(p.y, b->min.y, b->max.y);
    q.z = clampf(p.z, b->min.z, b->max.z);
    return q;
}

static int test(const Sphere* s, const Aabb* b, Contact* out) {
    if (s->r <= 0.0f) return -1;
    if (b->min.x > b->max.x || b->min.y > b->max.y || b->min.z > b->max.z) return -2;
    V3 q = closest(b, s->c);
    float dx = s->c.x - q.x, dy = s->c.y - q.y, dz = s->c.z - q.z;
    float d2 = dx*dx + dy*dy + dz*dz;
    out->hit = 0;
    out->depth = 0.0f;
    if (d2 > s->r * s->r) return 0;
    if (d2 > 1e-12f) {
        float d = sqrtf(d2);
        out->normal.x = dx / d; out->normal.y = dy / d; out->normal.z = dz / d;
        out->depth = s->r - d;
    } else {
        /* centre inside the box: push out along the shallowest axis */
        float px = fminf(s->c.x - b->min.x, b->max.x - s->c.x);
        float py = fminf(s->c.y - b->min.y, b->max.y - s->c.y);
        float pz = fminf(s->c.z - b->min.z, b->max.z - s->c.z);
        float m = px; int ax = 0;
        if (py < m) { m = py; ax = 1; }
        if (pz < m) { m = pz; ax = 2; }
        out->normal.x = out->normal.y = out->normal.z = 0.0f;
        if (ax == 0) out->normal.x = (s->c.x < (b->min.x + b->max.x) * 0.5f) ? -1.0f : 1.0f;
        if (ax == 1) out->normal.y = (s->c.y < (b->min.y + b->max.y) * 0.5f) ? -1.0f : 1.0f;
        if (ax == 2) out->normal.z = (s->c.z < (b->min.z + b->max.z) * 0.5f) ? -1.0f : 1.0f;
        out->depth = m + s->r;
    }
    out->point = q;
    out->hit = 1;
    return 1;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Aabb b = { {0,0,0}, {1,1,1} };
    Contact c;

    Sphere far = { {5,5,5}, 0.5f };
    check("miss_far", test(&far, &b, &c) == 0);

    Sphere touch = { {1.4f,0.5f,0.5f}, 0.5f };
    check("hit_face", test(&touch, &b, &c) == 1);
    check("normal_points_out", c.normal.x > 0.9f);
    check("depth_positive", c.depth > 0.0f && c.depth < 0.2f);

    /* exactly touching: distance == radius */
    Sphere edge = { {1.5f,0.5f,0.5f}, 0.5f };
    check("grazing_hits", test(&edge, &b, &c) == 1);

    Sphere inside = { {0.5f,0.5f,0.5f}, 0.25f };
    check("inside_hit", test(&inside, &b, &c) == 1);
    check("inside_pushes_out", c.depth > 0.0f);
    check("inside_normal_unit",
          fabsf(sqrtf(c.normal.x*c.normal.x + c.normal.y*c.normal.y + c.normal.z*c.normal.z) - 1.0f) < 1e-5f);

    Sphere corner = { {1.2f,1.2f,1.2f}, 0.5f };
    check("corner_hit", test(&corner, &b, &c) == 1);
    check("corner_normal_diagonal", c.normal.x > 0.5f && c.normal.y > 0.5f && c.normal.z > 0.5f);

    Sphere bad = { {0,0,0}, 0.0f };
    check("zero_radius_rejected", test(&bad, &b, &c) == -1);
    Aabb inv = { {1,1,1}, {0,0,0} };
    check("inverted_box_rejected", test(&touch, &inv, &c) == -2);

    /* Log line: depth=0.0031 normal=(0,1,0), radius=0.5, floor y=0. */
    Sphere logged_ball = { {0.0f, 0.4969f, 0.0f}, 0.5f };
    Contact logged_contact;
    memset(&logged_contact, 0, sizeof(logged_contact));
    logged_contact.normal = (V3){0.0f, 1.0f, 0.0f};
    logged_contact.depth = 0.0031f;
    DepthAudit good_depth = audit_plane_depth(
        logged_ball.c, (V3){0,0,0}, logged_ball.r, logged_contact, 5.1e-5f);
    check("logged_depth_self_consistent", good_depth.ok);
    check("logged_center_height_inverted",
          fabsf(good_depth.center_distance_from_depth - 0.4969f) < 1e-6f);

    Contact wrong_depth = logged_contact;
    wrong_depth.depth = 0.0042f;
    DepthAudit bad_depth = audit_plane_depth(
        logged_ball.c, (V3){0,0,0}, logged_ball.r, wrong_depth, 5.1e-5f);
    check("logged_depth_mismatch_detected", !bad_depth.ok && !bad_depth.depth_ok);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
