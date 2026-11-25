
/* roll: rolling resistance as a torque proportional to angular velocity.
 *
 * Model: tau = -c * omega, so I * dw/dt = -c * w and with rolling without
 * slipping (v = w R) the centre-of-mass speed decays exponentially:
 *
 *     v(t) = v0 * exp(-lambda t),   lambda = c / I
 *
 * This is NOT sliding friction: the contact point has no relative slip, the
 * loss comes from deformation, and the decay is exponential, not linear.
 *
 * Consequences used below:
 *   stop time to threshold v_th:  T = ln(v0 / v_th) / lambda
 *   stop distance to threshold:   D = (v0 - v_th) / lambda
 *   asymptotic distance (no threshold): D_inf = v0 / lambda
 *
 * Coefficient and threshold are solved backwards from design targets:
 *   lambda = ln(v_max / v_ref) / (t_max - t_target)
 *   v_th   = v_max * exp(-lambda * t_max)
 * so a ball at v_ref reads as stopped (below v_th) after t_target, and no
 * ball up to v_max can roll longer than t_max before the threshold snaps
 * it to zero. The snap exists because the exponential tail never reaches
 * zero on its own.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

typedef struct {
    float lambda;   /* decay rate c/I, 1/s */
    float v_th;     /* |v| below this is snapped to zero, m/s */
} Roll;

/* solid sphere I = (2/5) m R^2 */
static float sphere_inertia(float m, float r) { return 0.4f * m * r * r; }

/* torque coefficient c from decay rate and inertia */
static float roll_torque_coeff(const Roll* r, float inertia) {
    return r->lambda * inertia;
}

/* solve lambda and v_th from target/max stop times; see header comment */
static int roll_design(Roll* r, float v_ref, float t_target,
                       float v_max, float t_max) {
    if (v_ref <= 0.0f || v_max <= v_ref) return -1;
    if (t_target <= 0.0f || t_max <= t_target) return -1;
    r->lambda = logf(v_max / v_ref) / (t_max - t_target);
    r->v_th = v_max * expf(-r->lambda * t_max);
    return 0;
}

/* one exact step: v *= exp(-lambda dt), returns the (signed) displacement.
 * below the threshold the ball is snapped to rest and does not move. */
static float roll_step(const Roll* r, float* v, float dt) {
    float v0 = *v;
    if (dt < 0.0f) return 0.0f;
    if (fabsf(v0) <= r->v_th) { *v = 0.0f; return 0.0f; }
    float decay = expf(-r->lambda * dt);
    *v = v0 * decay;
    return v0 * (1.0f - decay) / r->lambda;   /* exact integral of v(t) */
}

/* analytic stop distance from v0 to the threshold */
static float roll_stop_distance(const Roll* r, float v0) {
    float a = fabsf(v0);
    if (a <= r->v_th) return 0.0f;
    return (a - r->v_th) / r->lambda;
}

/* analytic stop time from v0 to the threshold */
static float roll_stop_time(const Roll* r, float v0) {
    float a = fabsf(v0);
    if (a <= r->v_th) return 0.0f;
    return logf(a / r->v_th) / r->lambda;
}

/* simulate to rest; returns distance, writes stop time */
static float roll_simulate(const Roll* r, float v0, float dt, float* out_time) {
    float v = v0, x = 0.0f, t = 0.0f;
    while (v != 0.0f && t < 1e6f) {
        x += roll_step(r, &v, dt);
        t += dt;
    }
    if (out_time) *out_time = t;
    return fabsf(x);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* design: a 2 m/s ball should read as stopped in 2.5 s, and nothing
     * up to 10 m/s may roll longer than 4 s */
    Roll r;
    check("design_ok", roll_design(&r, 2.0f, 2.5f, 10.0f, 4.0f) == 0);
    printf("lambda=%.4f v_th=%.4f\n", r.lambda, r.v_th);

    float I = sphere_inertia(1.0f, 0.5f);
    printf("torque_coeff_c=%.5f (I=%.3f)\n", roll_torque_coeff(&r, I), I);

    /* trade-off sweep: stop distance grows like 1/lambda while the decay
     * time constant shrinks; big lambda stops crisply but looks sticky */
    printf("sweep: lambda tau=1/lambda D(v=10) T(v=10)\n");
    for (float k = 0.25f; k <= 4.01f; k *= 2.0f) {
        Roll s = { r.lambda * k, r.v_th };
        printf("sweep: %.4f %.3f %.3f %.3f\n",
               s.lambda, 1.0f / s.lambda,
               s.v_th < 10.0f ? 10.0f / s.lambda : 0.0f,
               roll_stop_time(&s, 10.0f));
    }

    /* three initial speeds: simulated stop distance vs analytic, <= 10% */
    float speeds[3] = { 1.5f, 4.0f, 9.0f };
    for (int i = 0; i < 3; i++) {
        float t = 0.0f;
        float d = roll_simulate(&r, speeds[i], 1.0f / 240.0f, &t);
        float da = roll_stop_distance(&r, speeds[i]);
        float dev = da > 0.0f ? fabsf(d - da) / da : 0.0f;
        char name[64];
        snprintf(name, sizeof(name), "stop_dist_v%.1f_within_10pct", speeds[i]);
        if (dev <= 0.10f) {
            check(name, 1);
        } else {
            printf("actual_distance_v%.1f=%.4f (analytic %.4f, dev %.1f%%)\n",
                   speeds[i], d, da, dev * 100.0f);
            check(name, 0);
        }
        printf("v0=%.1f sim_d=%.4f analytic_d=%.4f stop_t=%.3f\n",
               speeds[i], d, da, t);
    }

    /* max stop time bound holds at the design ceiling */
    float t_max_seen = 0.0f;
    roll_simulate(&r, 10.0f, 1.0f / 240.0f, &t_max_seen);
    check("stop_time_within_t_max", t_max_seen <= 4.0f + 1.0f / 240.0f);

    /* below the threshold the ball snaps to rest immediately */
    Roll snap = r;
    float v = r.v_th * 0.5f;
    float dx = roll_step(&snap, &v, 1.0f / 60.0f);
    check("below_threshold_snaps", v == 0.0f && dx == 0.0f);

    /* boundary: without the threshold the exponential tail never ends;
     * reaching 1e-3 m/s from 9 m/s would take far longer than t_max */
    float t_tail = logf(9.0f / 1e-3f) / r.lambda;
    printf("tail_time_to_1e-3=%.2f (t_max=4.0)\n", t_tail);
    check("tail_would_exceed_t_max", t_tail > 4.0f);

    Roll bad;
    check("bad_speeds", roll_design(&bad, 10.0f, 2.5f, 2.0f, 4.0f) == -1);
    check("bad_times", roll_design(&bad, 2.0f, 4.0f, 10.0f, 2.5f) == -1);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
