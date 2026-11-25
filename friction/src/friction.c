
/* friction: sphere-sphere impulse with Coulomb tangential clamp.
 *
 * Why the original pseudocode leaks on oblique impacts:
 *   1. no tangential impulse at all, so the slide that caused the oblique
 *      contact is untouched and the balls stay overlapped next frame;
 *   2. a fresh FULL normal impulse is recomputed every overlapping frame;
 *      while the balls slide around each other the contact normal rotates,
 *      part of the surviving tangential velocity reads as new approach
 *      velocity and the per-frame impulses run away (see `trace`);
 *   3. no approach guard and no positional projection.
 *
 * Solver here: the collision impulse is applied ONCE, on the first frame of
 * a contact episode (instantaneous-impulse / hard-disk model); later
 * overlapping frames only get a full positional projection along the
 * current normal, matching the correction style of the integrate/ module.
 *
 * Cap study:
 *   MODE_NAIVE   Jt = k * Jn applied blindly. Reversal appears exactly when
 *                k exceeds the stop-impulse ratio tan(theta)/(1+e).
 *   MODE_COULOMB Jt = min(stop-impulse, mu * Jn): the cap is derived from
 *                the impulse that removes the tangential velocity, so a
 *                sign reversal is impossible by construction.
 *
 * Sign convention: n points a -> b; a receives +(Jn n + Jt t), b the
 * opposite, t = unit tangent opposing the incoming slide.
 *
 * Usage:
 *   ./friction selftest
 *   ./friction sweep | baseline | final
 *   ./friction trace <scenario 0..4> <mode> <k> <guard> <poscorr> <legacy>
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define PI       3.14159265358979323846f
#define RADIUS   1.0f
#define INV_M    1.0f                  /* equal unit masses, K = inv_ma+inv_mb */
#define K_TOT    (2.0f * INV_M)
#define SPEED    10.0f
#define REST     0.5f
#define MU       0.5f
#define DT       0.02f
#define MAX_STEP 600
#define NSCEN    5
#define NK       41
#define REV_EPS  0.01f
#define SLOP     0.005f

typedef struct { float x, y, vx, vy; } Ball;
typedef struct { float x, y; } V2;

typedef enum { MODE_NAIVE = 0, MODE_COULOMB = 1 } Mode;

typedef struct {
    Mode  mode;
    float k;                          /* cap: |Jt| <= k * Jn */
    int   approach_guard;             /* 0 = original pseudocode behaviour */
    int   position_correct;
    int   legacy;                     /* 1 = re-apply a full impulse/frame */
} Cfg;

typedef struct {
    int   contact_frames;
    int   reversals;                  /* post-impulse slide flips sign */
    float slide;                      /* accumulated tangential travel */
    float max_pen;
    int   pass_through;               /* centre crossed the contact plane */
    float jn_first, jt_first;         /* first-frame impulses */
} Metrics;

static const float angles[NSCEN] = { 15.0f, 30.0f, 45.0f, 60.0f, 75.0f };
static int g_trace = 0;

static V2 sub2(float ax, float ay, float bx, float by) {
    V2 r = { ax - bx, ay - by };
    return r;
}
static float dot2(V2 a, V2 b) { return a.x * b.x + a.y * b.y; }
static float len2(V2 a) { return sqrtf(dot2(a, a)); }

/* b rests at origin; contact normal n0 = (1,0). a's velocity makes angle
 * theta with n0, so the incoming slide is SPEED*sin(theta) on t0=(0,1). */
static void setup(int s, Ball* a, Ball* b) {
    float rad = angles[s] * PI / 180.0f;
    memset(b, 0, sizeof(*b));
    a->vx = SPEED * cosf(rad);
    a->vy = SPEED * sinf(rad);
    float t_back = 6.0f * RADIUS / SPEED;
    a->x = -2.0f * RADIUS - a->vx * t_back;
    a->y = 0.0f            - a->vy * t_back;
}

static Metrics simulate(const Cfg* cfg, int s) {
    Metrics m;
    memset(&m, 0, sizeof(m));
    Ball a, b;
    setup(s, &a, &b);
    int active = 0;
    const float touch = 2.0f * RADIUS;

    for (int step = 0; step < MAX_STEP; step++) {
        Ball a0 = a, b0 = b;
        float u_hit = 0.0f;             /* time of contact within the frame */
        int rewound = 0;
        a.x += a.vx * DT;
        a.y += a.vy * DT;
        b.x += b.vx * DT;
        b.y += b.vy * DT;

        V2 d = sub2(b.x, b.y, a.x, a.y);
        float dist = len2(d);
        float pen = touch - dist;

        if (pen <= 0.0f) {
            active = 0;                 /* separation ends the episode */
            if (m.contact_frames > 0) break;
            continue;
        }

        /* continuous collision: on a new episode, rewind to the exact
         * contact instant |d(tau)| = 2r so the impulse uses the true
         * contact normal and pre-impact velocities. */
        if (!active && !cfg->legacy) {
            V2 d0 = sub2(b0.x, b0.y, a0.x, a0.y);
            V2 rv0 = sub2(b0.vx, b0.vy, a0.vx, a0.vy);
            float A = dot2(rv0, rv0);
            float Bv = 2.0f * dot2(d0, rv0);
            float C = dot2(d0, d0) - touch * touch;
            float disc = Bv * Bv - 4.0f * A * C;
            if (A > 1e-12f && disc > 0.0f) {
                float r1 = (-Bv - sqrtf(disc)) / (2.0f * A); /* time, seconds */
                if (r1 > 0.0f && r1 <= DT) {
                    u_hit = r1;
                    rewound = 1;
                }
            }
            if (rewound) {
                a.x = a0.x + a0.vx * u_hit;
                a.y = a0.y + a0.vy * u_hit;
                b.x = b0.x + b0.vx * u_hit;
                b.y = b0.y + b0.vy * u_hit;
                d = sub2(b.x, b.y, a.x, a.y);
                dist = len2(d);
                pen = touch - dist;
            }
        }

        V2 n;
        if (dist > 1e-9f) { n.x = d.x / dist; n.y = d.y / dist; }
        else { n.x = 1.0f; n.y = 0.0f; }

        V2 rv = sub2(b.vx, b.vy, a.vx, a.vy);
        float vn = dot2(rv, n);
        V2 tv = sub2(rv.x, rv.y, vn * n.x, vn * n.y);
        float tl = len2(tv);
        V2 t = { (tl > 1e-9f) ? tv.x / tl : 0.0f,
                 (tl > 1e-9f) ? tv.y / tl : 1.0f };

        /* impulse phase: first frame of an episode, or every frame in legacy */
        if (!active || cfg->legacy) {
            float jn = 0.0f, jt = 0.0f;
            if (vn < 0.0f || !cfg->approach_guard)
                jn = -(1.0f + REST) * vn / K_TOT;
            if (tl > 1e-9f && jn > 0.0f) {
                if (cfg->mode == MODE_NAIVE)
                    jt = cfg->k * jn;                /* blind fixed cap */
                else {
                    float stop = tl / K_TOT;          /* impulse killing slide */
                    jt = fminf(stop, MU * jn);        /* Coulomb-derived cap */
                }
            }
            a.vx += jn * n.x * INV_M + jt * t.x * INV_M;
            a.vy += jn * n.y * INV_M + jt * t.y * INV_M;
            b.vx -= jn * n.x * INV_M + jt * t.x * INV_M;
            b.vy -= jn * n.y * INV_M + jt * t.y * INV_M;

            if (!active) {
                m.jn_first = jn;
                m.jt_first = jt;
                active = 1;
            }
            /* reversal: post-impulse slide along the pre-impulse tangent */
            if (tl > REV_EPS) {
                V2 rv1 = sub2(b.vx, b.vy, a.vx, a.vy);
                if (rv1.x * t.x + rv1.y * t.y < -REV_EPS)
                    m.reversals++;
            }
            if (g_trace)
                fprintf(stderr,
                        "IMPULSE step=%d pen=%.4f n=(%.3f,%.3f) vn=%.4f "
                        "vt=%.4f jn=%.4f jt=%.4f\n",
                        step, pen, n.x, n.y, vn, dot2(tv, t), jn, jt);
        }

        /* after a CCD rewind, walk the remainder of the frame with the
         * post-impulse velocities so the step still spans the full DT */
        if (rewound) {
            float rem = DT - u_hit;
            a.x += a.vx * rem;
            a.y += a.vy * rem;
            b.x += b.vx * rem;
            b.y += b.vy * rem;
            d = sub2(b.x, b.y, a.x, a.y);
            dist = len2(d);
            pen = touch - dist;
            if (dist > 1e-9f) { n.x = d.x / dist; n.y = d.y / dist; }
        }

        /* positional phase: full projection of the overlap (hard disk) */
        if (cfg->position_correct && pen > SLOP) {
            float corr = pen - SLOP;
            a.x -= n.x * corr * 0.5f;
            a.y -= n.y * corr * 0.5f;
            b.x += n.x * corr * 0.5f;
            b.y += n.y * corr * 0.5f;
        }

        m.contact_frames++;
        m.slide += fabsf(dot2(tv, t)) * DT;
        if (pen > m.max_pen) m.max_pen = pen;
        /* tunnelling = centres approach closer than 1.5 radii; merely
         * deflecting around each other (d.x flips sign at full contact)
         * is legitimate hard-disk motion, not pass-through */
        if (dist < 1.5f * RADIUS) m.pass_through = 1;
    }
    return m;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(int argc, char** argv) {
    const char* cmd = argc > 1 ? argv[1] : "sweep";

    if (strcmp(cmd, "selftest") == 0) {
        Cfg coul = { MODE_COULOMB, MU, 1, 1, 0 };
        for (int s = 0; s < NSCEN; s++) {
            Metrics m = simulate(&coul, s);
            check("coulomb_no_reversal", m.reversals == 0);
            check("coulomb_no_passthrough", m.pass_through == 0);
            check("coulomb_cap_respected",
                  fabsf(m.jt_first) <= MU * m.jn_first + 1e-5f);
            /* first-contact penetration is bounded by one frame's travel;
             * the legacy re-impulse path instead grows past one radius */
            check("coulomb_pen_bounded", m.max_pen < 0.5f * RADIUS);
        }

        /* k=0 reproduces the pseudocode's split: no tangential impulse */
        Cfg none = { MODE_NAIVE, 0.0f, 1, 1, 0 };
        Metrics m0 = simulate(&none, 2);
        check("zero_friction_no_jt", fabsf(m0.jt_first) < 1e-6f);

        /* stop-impulse ratio for 60 deg is tan(60)/(1+e) ~= 1.15: below it
         * the slide merely slows, well above it the tangent velocity flips */
        Cfg under = { MODE_NAIVE, 1.00f, 1, 1, 0 };
        Cfg over  = { MODE_NAIVE, 1.50f, 1, 1, 0 };
        check("under_stop_no_reversal", simulate(&under, 3).reversals == 0);
        check("over_stop_reverses", simulate(&over, 3).reversals == 1);

        /* legacy per-frame full impulse on the grazing case goes unstable */
        Cfg legacy = { MODE_NAIVE, 0.0f, 1, 0, 1 };
        Metrics ml = simulate(&legacy, 4);
        check("legacy_reimpulse_penetrates", ml.max_pen > 0.5f);

        printf("all_ok=%d\n", failures == 0);
        return failures == 0 ? 0 : 1;
    }

    if (strcmp(cmd, "trace") == 0) {
        int s = argc > 2 ? atoi(argv[2]) : 4;
        Mode mode = argc > 3 ? (Mode)atoi(argv[3]) : MODE_COULOMB;
        float k = argc > 4 ? atof(argv[4]) : MU;
        Cfg cfg = { mode, k, argc > 5 ? atoi(argv[5]) : 1,
                    argc > 6 ? atoi(argv[6]) : 1,
                    argc > 7 ? atoi(argv[7]) : 0 };
        g_trace = 1;
        Metrics m = simulate(&cfg, s);
        fprintf(stderr, "RESULT frames=%d rev=%d slide=%.4f maxpen=%.4f "
                "through=%d jn=%.4f jt=%.4f ratio=%.4f\n",
                m.contact_frames, m.reversals, m.slide, m.max_pen,
                m.pass_through, m.jn_first, m.jt_first,
                m.jn_first > 0.0f ? m.jt_first / m.jn_first : 0.0f);
        return 0;
    }

    if (strcmp(cmd, "baseline") == 0) {
        printf("config,scenario_deg,contact_frames,reversals,slide,"
               "max_pen,pass_through\n");
        for (int s = 0; s < NSCEN; s++) {
            Cfg bug = { MODE_NAIVE, 0.0f, 0, 0, 1 };      /* original */
            Cfg guarded = { MODE_NAIVE, 0.0f, 1, 0, 1 };
            Metrics mb = simulate(&bug, s);
            Metrics mg = simulate(&guarded, s);
            printf("original_pseudocode,%g,%d,%d,%.6f,%.6f,%d\n",
                   angles[s], mb.contact_frames, mb.reversals, mb.slide,
                   mb.max_pen, mb.pass_through);
            printf("approach_guard_only,%g,%d,%d,%.6f,%.6f,%d\n",
                   angles[s], mg.contact_frames, mg.reversals, mg.slide,
                   mg.max_pen, mg.pass_through);
        }
        return 0;
    }

    if (strcmp(cmd, "final") == 0) {
        printf("scenario_deg,contact_frames,reversals,slide,max_pen,"
               "pass_through,jn_first,jt_first,jt_over_jn\n");
        for (int s = 0; s < NSCEN; s++) {
            Cfg cfg = { MODE_COULOMB, MU, 1, 1, 0 };
            Metrics m = simulate(&cfg, s);
            printf("%g,%d,%d,%.6f,%.6f,%d,%.6f,%.6f,%.6f\n",
                   angles[s], m.contact_frames, m.reversals, m.slide,
                   m.max_pen, m.pass_through, m.jn_first, m.jt_first,
                   m.jt_first / m.jn_first);
        }
        return 0;
    }

    /* sweep: fixed-cap grid over all five scenarios + the Coulomb point */
    printf("mode,k,scenario_deg,contact_frames,reversals,slide,max_pen,"
           "pass_through,jn_first,jt_first,jt_over_jn\n");
    for (int i = 0; i < NK; i++) {
        float k = 0.1f * i;                           /* 0.0 .. 4.0 */
        for (int s = 0; s < NSCEN; s++) {
            Cfg cfg = { MODE_NAIVE, k, 1, 1, 0 };
            Metrics m = simulate(&cfg, s);
            printf("naive,%.3f,%g,%d,%d,%.6f,%.6f,%d,%.6f,%.6f,%.6f\n",
                   k, angles[s], m.contact_frames, m.reversals, m.slide,
                   m.max_pen, m.pass_through, m.jn_first, m.jt_first,
                   m.jt_first / m.jn_first);
        }
    }
    for (int s = 0; s < NSCEN; s++) {
        Cfg cfg = { MODE_COULOMB, MU, 1, 1, 0 };
        Metrics m = simulate(&cfg, s);
        printf("coulomb,%.3f,%g,%d,%d,%.6f,%.6f,%d,%.6f,%.6f,%.6f\n",
               MU, angles[s], m.contact_frames, m.reversals, m.slide,
               m.max_pen, m.pass_through, m.jn_first, m.jt_first,
               m.jt_first / m.jn_first);
    }
    return 0;
}
