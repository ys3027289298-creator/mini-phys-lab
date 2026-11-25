
/* rollfeel: translate the art note "rolls like ice" into tunable metrics.
 *
 * Model/test protocol (all explicit because feel bands are calibration data):
 *   - solid sphere, I = kappa*m*R^2, kappa = 2/5;
 *   - dimensionless rolling-resistance coefficient cr, F_r = cr*m*g;
 *   - Coulomb friction coefficient mu;
 *   - standardized slightly under-rotated start u0 = eps*v0;
 *   - rate/slip are sampled after h seconds;
 *   - energy ratio is measured over one reference metre at speed v*.
 *
 * With no arguments the program runs its self-checks. With two arguments it
 * evaluates one (mu, cr) pair and prints the actual values plus violation
 * directions.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define KAPPA 0.4f

typedef struct {
    float g;       /* gravity, m/s^2 */
    float v0;      /* reference/initial centre speed, m/s */
    float eps;     /* initial contact slip as a fraction of v0 */
    float h;       /* short measurement window, s */
    float length;  /* energy-ratio reference distance, m */
} RollConfig;

typedef struct {
    float d_rate;      /* D_h: fractional speed-loss rate over h, s^-1 */
    float slip;        /* S_h: contact slip speed / centre speed at h */
    float work_ratio;  /* eta_L: rolling-resistance work / rolling KE */
    float grip_time;   /* time to reach zero slip, s; INFINITY if never */
} Metrics;

typedef struct {
    float d_lo, d_hi;
    float s_lo, s_hi;
    float e_lo, e_hi;
} Targets;

static const RollConfig DEFAULT_CFG = {
    9.81f, 5.0f, 0.10f, 0.020f, 1.0f
};

/* Art calibration used here as a concrete, replaceable default. */
static const Targets DEFAULT_TARGETS = {
    0.55f, 1.15f,   /* D_h, s^-1 */
    0.020f, 0.050f, /* S_h */
    0.020f, 0.060f  /* eta_1m */
};

static int evaluate(float mu, float cr, const RollConfig* cfg, Metrics* out) {
    if (!isfinite(mu) || !isfinite(cr)) return -1;
    if (mu < 0.0f || mu > 1.0f || cr < 0.0f || cr > 1.0f) return -1;
    if (cfg->g <= 0.0f || cfg->v0 <= 0.0f || cfg->h <= 0.0f || cfg->length <= 0.0f)
        return -2;
    if (cfg->eps < 0.0f || cfg->eps >= 1.0f) return -2;

    const float A = 1.0f + 1.0f / KAPPA; /* coefficient on mu in slip closure */
    const float u0 = cfg->eps * cfg->v0;
    const float closure = cfg->g * (A * mu + cr);

    out->grip_time = closure > 1e-12f ? u0 / closure : INFINITY;

    if (out->grip_time >= cfg->h) {
        /* Still sliding: vdot = -g(mu+cr), udot = -g(A*mu+cr). */
        out->d_rate = cfg->g * (mu + cr) / cfg->v0;

        float uh = u0 - cfg->g * (A * mu + cr) * cfg->h;
        float vh = cfg->v0 - cfg->g * (mu + cr) * cfg->h;
        out->slip = uh > 0.0f ? uh / fmaxf(vh, 1e-6f) : 0.0f;
    } else {
        /* Friction impulse has already established rolling. Rolling
         * deceleration is cr*g/(1+kappa) after the static-friction torque
         * is included. */
        float ar = cr * cfg->g / (1.0f + KAPPA);
        float sliding_loss = cfg->g * (mu + cr) * out->grip_time;
        float rolling_loss = ar * (cfg->h - out->grip_time);
        out->d_rate = (sliding_loss + rolling_loss) / (cfg->v0 * cfg->h);
        out->slip = 0.0f;
    }

    /* W_r = cr*m*g*L; K = 1/2*m*(1+kappa)*v0^2. */
    out->work_ratio =
        2.0f * cr * cfg->g * cfg->length /
        ((1.0f + KAPPA) * cfg->v0 * cfg->v0);
    return 0;
}

static int in_range(float x, float lo, float hi) {
    return isfinite(x) && x >= lo && x <= hi;
}

static int all_targets_met(const Metrics* m, const Targets* t) {
    return in_range(m->d_rate, t->d_lo, t->d_hi) &&
           in_range(m->slip, t->s_lo, t->s_hi) &&
           in_range(m->work_ratio, t->e_lo, t->e_hi);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static void report_metric(const char* name, float value, float lo, float hi,
                          const char* low_hint, const char* high_hint) {
    if (value < lo) {
        printf("%s=%.6g actual<[%.6g,%.6g] direction=low %s\n",
               name, value, lo, hi, low_hint);
    } else if (value > hi) {
        printf("%s=%.6g actual>[%.6g,%.6g] direction=high %s\n",
               name, value, lo, hi, high_hint);
    } else {
        printf("%s=%.6g in_[%.6g,%.6g] status=ok\n", name, value, lo, hi);
    }
}

static int report_pair(float mu, float cr) {
    Metrics m;
    int rc = evaluate(mu, cr, &DEFAULT_CFG, &m);
    if (rc != 0) {
        printf("evaluate failed rc=%d for mu=%g cr=%g\n", rc, mu, cr);
        return 2;
    }

    printf("mu=%.4f cr=%.4f grip_time=%.4fs\n", mu, cr,
           isinf(m.grip_time) ? -1.0f : m.grip_time);
    report_metric("D_20ms", m.d_rate, DEFAULT_TARGETS.d_lo, DEFAULT_TARGETS.d_hi,
                  "(减速太慢, 像冰)", "(减速太突兀, 摩擦/毡感过重)");
    report_metric("S_20ms", m.slip, DEFAULT_TARGETS.s_lo, DEFAULT_TARGETS.s_hi,
                  "(低于下限: 过于玻璃/冰感, 或被冲量瞬间抓死)",
                  "(残余滑移过高, 抓不住/像滑面)");
    report_metric("eta_1m", m.work_ratio, DEFAULT_TARGETS.e_lo, DEFAULT_TARGETS.e_hi,
                  "(阻力功占比过低, 滚动不费力/像冰)",
                  "(阻力功占比过高, 像沙/很快停)");

    int ok = all_targets_met(&m, &DEFAULT_TARGETS);
    printf("all_metrics_in_target=%d\n", ok);
    return ok ? 0 : 1;
}

int main(int argc, char** argv) {
    if (argc == 3) {
        char* end_mu = NULL;
        char* end_cr = NULL;
        float mu = strtof(argv[1], &end_mu);
        float cr = strtof(argv[2], &end_cr);
        if (*end_mu != '\0' || *end_cr != '\0') {
            fprintf(stderr, "usage: %s [mu cr]\n", argv[0]);
            return 2;
        }
        return report_pair(mu, cr);
    }

    Metrics m;

    /* Recommended interior point and all four corners of the conservative
     * all-pairs tuning box mu in [0.36,0.49], cr in [0.05,0.09]. */
    check("recommended_eval_ok", evaluate(0.40f, 0.07f, &DEFAULT_CFG, &m) == 0);
    check("recommended_in_target", all_targets_met(&m, &DEFAULT_TARGETS));
    check("recommended_d_value", fabsf(m.d_rate - 0.92214f) < 1e-4f);
    check("recommended_s_value", fabsf(m.slip - 0.0431123f) < 1e-4f);
    check("recommended_eta_value", fabsf(m.work_ratio - 0.03924f) < 1e-5f);

    static const float box_mu[4] = {0.36f, 0.49f, 0.36f, 0.49f};
    static const float box_cr[4] = {0.05f, 0.05f, 0.09f, 0.09f};
    for (int i = 0; i < 4; ++i) {
        Metrics corner;
        char name[64];
        snprintf(name, sizeof(name), "safe_box_corner_%d", i);
        check(name, evaluate(box_mu[i], box_cr[i], &DEFAULT_CFG, &corner) == 0 &&
                    all_targets_met(&corner, &DEFAULT_TARGETS));
    }

    /* Ice pair: slow decay, too much residual slip after the standardized
     * under-rotation, and too little rolling-resistance work. */
    Metrics ice;
    evaluate(0.03f, 0.008f, &DEFAULT_CFG, &ice);
    check("ice_d_too_low", ice.d_rate < DEFAULT_TARGETS.d_lo);
    check("ice_s_too_high_in_disturbance_test", ice.slip > DEFAULT_TARGETS.s_hi);
    check("ice_eta_too_low", ice.work_ratio < DEFAULT_TARGETS.e_lo);

    /* Over-grippy/high-resistance pair. */
    Metrics harsh;
    evaluate(0.9f, 0.15f, &DEFAULT_CFG, &harsh);
    check("harsh_d_too_high", harsh.d_rate > DEFAULT_TARGETS.d_hi);
    check("harsh_s_too_low", harsh.slip < DEFAULT_TARGETS.s_lo);
    check("harsh_eta_too_high", harsh.work_ratio > DEFAULT_TARGETS.e_hi);

    check("bad_mu_rejected", evaluate(1.5f, 0.1f, &DEFAULT_CFG, &m) == -1);
    check("bad_cr_rejected", evaluate(0.4f, -0.1f, &DEFAULT_CFG, &m) == -1);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
