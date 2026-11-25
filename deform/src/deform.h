/*
 * deform: compliant sphere contact (squash & bounce).
 *
 * Model
 * -----
 * During contact the sphere is squeezed against a rigid floor by a
 * spring-damper along the contact normal.  Let x >= 0 be the compression
 * (the amount the vertical diameter is shortened, in metres) and v_impact
 * the centre speed at first touch, directed into the floor.
 *
 *   m x'' = -k x - c x'          (x > 0, contact; free flight otherwise)
 *
 * Material stiffness k [N/m] and damping c [N s/m].
 *
 * Closed-form answers used for calibration (under-damped contact, zeta < 1):
 *
 *   natural frequency   wn   = sqrt(k / m)
 *   damping ratio       zeta = c / (2 m wn)
 *   damped frequency    wd   = wn sqrt(1 - zeta^2)
 *
 * Energy balance at maximum compression (velocity x' = 0) gives the
 * squash amplitude as a function of impact speed and stiffness:
 *
 *   x_max = A(zeta) * v_impact / wn,
 *
 *   A(zeta) = exp(-zeta * atan(sqrt(1-zeta^2)/zeta) / sqrt(1-zeta^2))
 *
 * so with zeta fixed, x_max = A v_impact sqrt(m / k).
 *
 * Max squash ratio (the quantity measured frame-by-frame off the video,
 * vertical shortening divided by the resting diameter 2r):
 *
 *   epsilon = x_max / (2r)
 *
 * Invert it to read stiffness off a target squash ratio epsilon_target:
 *
 *   k = m * (A v_impact / (2 r epsilon_target))^2
 *
 * Bounce height: the same damped oscillator leaves contact at speed
 * v_out = e v_impact with coefficient of restitution
 *
 *   e = exp(-pi zeta / sqrt(1 - zeta^2))
 *
 * hence h_bounce = e^2 h_drop.  Holding zeta (i.e. e) fixed while sweeping
 * k keeps the bounce-height curve flat: stiffness controls squash
 * visibility, restitution controls bounce, and the two stay decoupled.
 */
#ifndef DEFORM_H
#define DEFORM_H

#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    float mass;        /* kg           */
    float radius;      /* m            */
    float stiffness;   /* k, N/m       */
    float damping;     /* c, N s/m     */
    float gravity;     /* m/s^2        */
} DeformBody;

typedef struct {
    float impact_speed;   /* centre speed at first touch, m/s  */
    float x_max;          /* peak compression, m                */
    float squash_ratio;   /* x_max / (2r)                       */
    float exit_speed;     /* centre speed leaving contact, m/s  */
    float bounce_height;  /* centre rise after contact, m       */
    float contact_time;   /* s, 0 if contact never resolved     */
} ContactResult;

/* ---- parameters ------------------------------------------------------ */

static inline int deform_valid(const DeformBody* b) {
    return b != 0 &&
           b->mass > 0.0f && b->radius > 0.0f &&
           b->stiffness > 0.0f && b->damping >= 0.0f &&
           b->gravity > 0.0f;
}

static inline float deform_wn(const DeformBody* b) {
    return sqrtf(b->stiffness / b->mass);
}

static inline float deform_zeta(const DeformBody* b) {
    float wn = deform_wn(b);
    return b->damping / (2.0f * b->mass * wn);
}

/* coefficient of restitution of the spring-damper contact */
static inline float deform_restitution(float zeta) {
    if (zeta >= 1.0f) return 0.0f;
    return expf(-(float)M_PI * zeta / sqrtf(1.0f - zeta * zeta));
}

/* damping ratio that produces the requested coefficient of restitution */
static inline float deform_zeta_from_restitution(float e) {
    if (e <= 0.0f) return 0.999999f;
    if (e >= 1.0f) return 0.0f;
    float l = -logf(e);
    return l / sqrtf((float)M_PI * (float)M_PI + l * l);
}

/*
 * Peak-compression factor A(zeta): x_max = A(zeta) * v_impact / wn.
 * Undamped limit A(0) = 1, i.e. pure energy balance k x^2/2 = m v^2/2.
 */
static inline float deform_peak_factor(float zeta) {
    if (zeta <= 0.0f) return 1.0f;
    if (zeta >= 1.0f) return 0.0f;
    float s = sqrtf(1.0f - zeta * zeta);
    return expf(-zeta * atanf(s / zeta) / s);
}

/*
 * Squash amplitude from impact speed and material stiffness.
 * Returns -1 on bad inputs or an over-damped (zeta >= 1) configuration,
 * where the closed-form under-damped peak does not apply.
 */
static inline float deform_x_max(float v_impact, const DeformBody* b) {
    if (!deform_valid(b) || v_impact <= 0.0f) return -1.0f;
    float zeta = deform_zeta(b);
    if (zeta >= 1.0f) return -1.0f;
    return deform_peak_factor(zeta) * v_impact / deform_wn(b);
}

/*
 * Inverse calibration: stiffness [N/m] that produces target squash ratio
 * epsilon_target at the given impact speed, holding zeta fixed.
 * Returns -1 on bad inputs.
 */
static inline float deform_stiffness_for(float v_impact, float epsilon_target,
                                         float mass, float radius, float zeta) {
    if (v_impact <= 0.0f || epsilon_target <= 0.0f ||
        mass <= 0.0f || radius <= 0.0f || zeta < 0.0f || zeta >= 1.0f)
        return -1.0f;
    float a = deform_peak_factor(zeta);
    float wn = a * v_impact / (2.0f * radius * epsilon_target);
    return mass * wn * wn;
}

/* damping coefficient c for a target restitution e at stiffness k */
static inline float deform_damping_for(float mass, float stiffness, float e) {
    float zeta = deform_zeta_from_restitution(e);
    return 2.0f * mass * zeta * sqrtf(stiffness / mass);
}

/*
 * Visual squash & stretch scales with volume (area) preservation:
 * vertical scale sy, horizontal scale sx = 1/sqrt(sy).
 * sy in (0, 1] during compression; sx >= 1.
 */
static inline void deform_scales(float squash_ratio, float* sx, float* sy) {
    float s = 1.0f - fmaxf(0.0f, squash_ratio);
    if (s <= 0.0f) s = 1e-4f;
    if (sy) *sy = s;
    if (sx) *sx = 1.0f / sqrtf(s);
}

/*
 * Time-stepped contact simulation (semi-implicit Euler, matching the
 * integrate module).  The body starts at first touch: centre height r,
 * velocity -v_impact, and runs until it leaves contact (centre back at
 * r with upward velocity).  Compression x = max(0, r - centre_height);
 * the spring law acts only while x > 0.
 *
 * Bounce height follows from the exit speed ballistically: v^2 / (2g).
 * dt must resolve the contact oscillation: dt << 1/wn.  Returns -1 on
 * bad inputs; contact_time == 0 in the result means contact never
 * resolved within max_steps (e.g. over-damped).
 */
static inline int deform_simulate(const DeformBody* b, float v_impact,
                                  float dt, int max_steps, ContactResult* out) {
    if (!deform_valid(b) || v_impact <= 0.0f || dt <= 0.0f ||
        max_steps <= 0 || out == 0)
        return -1;

    float y = b->radius;          /* centre height, floor at y = 0 */
    float vy = -v_impact;
    float x_max = 0.0f;
    float contact_time = 0.0f;
    int contacted = 0, resolved = 0;

    for (int i = 0; i < max_steps; i++) {
        float x = fmaxf(0.0f, b->radius - y);
        /* compression rate x' = -vy; contact ends once the centre
         * rebounds past the resting height with upward velocity */
        if (contacted && y >= b->radius && vy > 0.0f) { resolved = 1; break; }

        float xdot = -vy;
        /* y-frame forces: spring and damper both push upward while the
         * ball is compressing (xdot > 0), damper reverses during rebound */
        float f_spring = b->stiffness * x;
        float f_damp = (x > 0.0f) ? b->damping * xdot : 0.0f;
        float acc = -b->gravity + (f_spring + f_damp) / b->mass;
        vy += acc * dt;
        y += vy * dt;
        contact_time += dt;

        if (x > 0.0f) {
            contacted = 1;
            if (x > x_max) x_max = x;
        }
    }

    out->impact_speed = v_impact;
    out->x_max = x_max;
    out->squash_ratio = x_max / (2.0f * b->radius);
    out->exit_speed = resolved ? fmaxf(0.0f, vy) : 0.0f;
    out->bounce_height = resolved ? vy * vy / (2.0f * b->gravity) : 0.0f;
    out->contact_time = resolved ? contact_time : 0.0f;
    return 0;
}

#endif /* DEFORM_H */
