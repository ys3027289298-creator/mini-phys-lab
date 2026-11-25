/*
 * degenerate self-tests.
 *
 * Every engineered exact-coincidence case must come out with a finite unit
 * normal; failures name the ball pair explicitly (ball A id / ball B id).
 */
#include "degenerate.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

static int failures = 0;
static void check(const char *name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* One exact-coincidence case; returns 0 only if the normal is unusable. */
static int coincidence_case(const char *tag, int id_a, int id_b,
                            const Sphere *a, const Sphere *b,
                            V3 rel_vel, PairHistory *h,
                            NormalSource expect) {
    Contact c = sphere_contact(a, b, rel_vel, h);
    int ok = c.hit && v3_is_unit(c.normal, 1e-5f) && c.source == expect;
    if (!ok) {
        printf("  FAIL pair (ball %d, ball %d) tag=%s: source=%d "
               "normal=(%.6g,%.6g,%.6g) bad=%d unit=%d\n",
               id_a, id_b, tag, c.source,
               c.normal.x, c.normal.y, c.normal.z,
               v3_is_bad(c.normal), v3_is_unit(c.normal, 1e-5f));
    }
    return ok;
}

int main(void) {
    /* ---- ordinary separated / overlapping pairs ---------------------- */
    Sphere a = {{0, 0, 0}, 1.0f};
    Sphere b = {{1.5f, 0, 0}, 1.0f};
    Contact c = sphere_contact(&a, &b, V3_ZERO, NULL);
    check("overlap_hit", c.hit == 1);
    check("position_normal", c.source == NORMAL_FROM_POSITION
          && fabsf(c.normal.x + 1.0f) < 1e-6f);
    check("overlap_depth", fabsf(c.depth - 0.5f) < 1e-6f);

    Sphere far = {{5, 0, 0}, 1.0f};
    c = sphere_contact(&a, &far, V3_ZERO, NULL);
    check("separated_no_hit", c.hit == 0);

    /* ---- coincidence level 2: velocity resolves the limit ------------- */
    Sphere same1 = {{2, 2, 2}, 1.0f};
    Sphere same2 = {{2, 2, 2}, 1.0f};
    PairFrame storage[8];
    PairHistory h; history_init(&h, storage, 8);
    V3 head_on = {-3, 4, 0};           /* B moving toward A, nonzero */
    check("coinc_velocity",
          coincidence_case("vel", 1, 2, &same1, &same2, head_on, &h,
                           NORMAL_FROM_VELOCITY));
    c = sphere_contact(&same1, &same2, head_on, &h);
    check("velocity_normal_heading",
          fabsf(c.normal.x - 0.6f) < 1e-5f && fabsf(c.normal.y + 0.8f) < 1e-5f);

    /* ---- coincidence level 3: dead still now, history carries motion -- */
    V3 still = {0, 0, 0};
    history_push(&h, (V3){0.10f, 0, 0}, (V3){1, 0, 0});   /* age 1 */
    history_push(&h, (V3){0.05f, 0, 0}, (V3){1, 0, 0});   /* age 2 */
    check("coinc_history_position",
          coincidence_case("hist_pos", 1, 2, &same1, &same2, still, &h,
                           NORMAL_FROM_HISTORY));
    c = sphere_contact(&same1, &same2, still, &h);
    check("history_normal_points_past", c.normal.x < -0.99f);
    check("history_age_reported", c.history_age >= 1 && c.history_age <= 2);

    /* velocity-only history: positions stayed coincident, motion recorded */
    PairFrame vstore[4]; PairHistory vh; history_init(&vh, vstore, 4);
    history_push(&vh, V3_ZERO, (V3){0, -2, 0});
    check("coinc_history_velocity",
          coincidence_case("hist_vel", 3, 4, &same1, &same2, still, &vh,
                           NORMAL_FROM_HISTORY));
    c = sphere_contact(&same1, &same2, still, &vh);
    check("history_vel_normal", fabsf(c.normal.y - 1.0f) < 1e-5f);

    /* empty history window: fixed axis is the last resort, never zero */
    PairFrame estore[2]; PairHistory eh; history_init(&eh, estore, 2);
    check("coinc_fixed_axis",
          coincidence_case("fixed", 5, 6, &same1, &same2, still, &eh,
                           NORMAL_FROM_FIXED_AXIS));
    /* NULL history means "no buffer configured": same fixed-axis outcome */
    check("no_history_same_fixed",
          coincidence_case("fixed_null", 5, 6, &same1, &same2, still, NULL,
                           NORMAL_FROM_FIXED_AXIS));
    /* ---- naive formula crashes at coincidence (the original bug) ------ */
    Contact bad = sphere_contact_naive(&same1, &same2);
    check("naive_normal_is_nan", v3_is_bad(bad.normal));
    V3 va = {0, 0, 0}, vb = {0, 0, 0};
    apply_impulse(&va, &vb, bad.normal, 0.5f);
    check("naive_impulse_is_nan", v3_is_bad(va) && v3_is_bad(vb));

    /* ---- fixed-axis response slides, never bounces -------------------- */
    /*
     * Two coincident spheres with relative motion orthogonal to the fixed
     * axis. The impulse along x is then zero: no kinetic energy is removed
     * and nothing bounces -- the pair only keeps sliding in y, drifting apart
     * along the axis only once positional correction moves the centres.
     */
    Sphere s1 = {{0, 0, 0}, 1.0f}, s2 = {{0, 0, 0}, 1.0f};
    V3 axis = {1, 0, 0};

    /* FIXED_AXIS needs the whole chain exhausted: zero position difference,
     * zero relative velocity, and an empty history window. */
    PairFrame fstore[2]; PairHistory fh; history_init(&fh, fstore, 2);
    Contact fc = sphere_contact(&s1, &s2, V3_ZERO, &fh);
    check("fixed_axis_is_x", fc.source == NORMAL_FROM_FIXED_AXIS
          && fabsf(fc.normal.x - 1.0f) < 1e-6f);

    /* Motion is purely tangential to the chosen axis: impulse is a no-op. */
    V3 v1 = {0, 0, 0}, v2 = {0, 5, 0};
    float ke_before = 0.5f * (v1.x*v1.x + v1.y*v1.y + v1.z*v1.z
                            + v2.x*v2.x + v2.y*v2.y + v2.z*v2.z);
    apply_impulse(&v1, &v2, axis, 0.5f);
    float ke_after = 0.5f * (v1.x*v1.x + v1.y*v1.y + v1.z*v1.z
                           + v2.x*v2.x + v2.y*v2.y + v2.z*v2.z);
    check("fixed_axis_no_energy_loss", fabsf(ke_after - ke_before) < 1e-5f);
    check("fixed_axis_slides_y", fabsf(v2.y - 5.0f) < 1e-5f
          && fabsf(v1.y) < 1e-5f);
    check("fixed_axis_no_bounce_x", fabsf(v2.x) < 1e-5f
          && fabsf(v1.x) < 1e-5f);

    /* contrast: a velocity-resolved normal turns the same motion into a
     * genuine bounce (energy redirected along the contact normal) */
    V3 w1 = {0, 0, 0}, w2 = {0, 5, 0};
    Contact gc = sphere_contact(&s1, &s2, (V3){0, 5, 0}, &fh);
    apply_impulse(&w1, &w2, gc.normal, 0.5f);
    check("velocity_path_bounces", gc.source == NORMAL_FROM_VELOCITY
          && w1.y > 0.0f && fabsf(w1.y - 3.75f) < 1e-4f
          && fabsf(w2.y - 1.25f) < 1e-4f);

    /* ---- ring buffer wraps and reports oldest-first correctly --------- */
    PairFrame rstore[3]; PairHistory rh; history_init(&rh, rstore, 3);
    for (int k = 0; k < 10; k++)
        history_push(&rh, (V3){(float)k, 0, 0}, V3_ZERO);
    check("ring_saturates", rh.count == 3);
    /* newest pushed was k=9 at +x, so the B->A normal must be -x */
    Contact rc = sphere_contact(&s1, &s2, V3_ZERO, &rh);
    check("ring_newest_first", rc.source == NORMAL_FROM_HISTORY
          && rc.normal.x < -0.99f && rc.history_age == 1);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
