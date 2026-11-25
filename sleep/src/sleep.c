
/* sleep: velocity sleep policy for settled balls.
 *
 * Criterion (both components must pass for confirm_frames consecutive steps):
 *     horizontal speed sqrt(vx^2 + vz^2) < v_horizontal
 *     vertical speed   |vy|             < v_vertical
 * Checking both is mandatory: a ball with vy == 0 that is still travelling
 * horizontally must never be put to sleep.
 *
 * Threshold selection experiment:
 *   - isolated rolling sentinels on each of three grades give the false-sleep
 *     curve and validate the analytic minimum rolling speed;
 *   - a flat settling pile gives the sleep-ratio / collision-skip curve;
 *   - a never-sleep oracle audits every sleep event by future displacement.
 */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* ============================== math types ============================== */

typedef struct { float x, y, z; } V3;

static float vdot(V3 a, V3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static V3    vsub(V3 a, V3 b) { V3 r = { a.x-b.x, a.y-b.y, a.z-b.z }; return r; }
static V3    vadd(V3 a, V3 b) { V3 r = { a.x+b.x, a.y+b.y, a.z+b.z }; return r; }
static V3    vscale(V3 a, float s) { V3 r = { a.x*s, a.y*s, a.z*s }; return r; }

/* ============================== sleep policy =========================== */

typedef struct {
    float v_horizontal;   /* max horizontal speed, m/s */
    float v_vertical;     /* max |vertical speed|,  m/s */
    int   confirm_frames; /* consecutive slow frames required */
} SleepCfg;

typedef struct {
    int slow_count;
    int asleep;
} SleepState;

/* Component-wise gate. Horizontal plane is xz; vertical axis is y. */
static int sleep_is_slow(const V3* v, const SleepCfg* c) {
    float vh = sqrtf(v->x*v->x + v->z*v->z);
    float vv = fabsf(v->y);
    return vh < c->v_horizontal && vv < c->v_vertical;
}

/* One evaluation per frame, called AFTER integration + solving.
 * Returns 1 on the awake -> asleep transition. */
static int sleep_update(SleepState* s, const V3* v, const SleepCfg* c) {
    if (s->asleep) return 0;
    if (sleep_is_slow(v, c)) {
        s->slow_count++;
        if (s->slow_count >= c->confirm_frames) {
            s->asleep = 1;
            return 1;
        }
    } else {
        s->slow_count = 0;
    }
    return 0;
}

static void sleep_wake(SleepState* s) {
    s->asleep = 0;
    s->slow_count = 0;
}

static void sleep_init(SleepState* s) { memset(s, 0, sizeof(*s)); }

/* ====================== slope derivation constants ===================== */

#define G_MS2            9.8f
#define DT_S             (1.0f / 60.0f)
#define CONFIRM_FRAMES   30            /* 0.5 s confirmation window */
#define SOLID_SPHERE_K   (5.0f / 7.0f) /* rolling acceleration factor */
#define ROLL_RESIST      0.03f         /* mu_r, must be < tan(theta_min) */
#define LIN_DAMP_PER_S   0.5f          /* mild viscous/air damping */

/* Three slope grades, degrees. */
#define GRADE_COUNT 3
static const float GRADES_DEG[GRADE_COUNT] = { 5.0f, 10.0f, 20.0f };

/* Operating threshold picked from the sweep (artifacts/sleep_sweep.csv):
 * above the solver normal-jitter floor (~0.016) and ~3x below the shallow
 * grade confirmation-window minimum rolling speed (~0.14). */
#define THRESHOLD_SELECTED 0.045f

#define PI_F 3.14159265358979323846f

/* Net downhill acceleration of a uniform solid sphere on a grade.
 * Gravity's tangential kick is partly converted to spin (factor 5/7),
 * while rolling resistance opposes translation at full mu_r*g*cos(theta):
 *   a = (5/7) g sin theta - mu_r g cos theta. */
static float slope_accel(float deg) {
    float th = deg * PI_F / 180.0f;
    return SOLID_SPHERE_K * G_MS2 * sinf(th)
           - ROLL_RESIST * G_MS2 * cosf(th);
}

/* Speed on the Nth evaluated frame from rest, using the exact discrete
 * recurrence the simulator applies (damping last):
 *     v_{k+1} = q (v_k + a dt),  q = exp(-c dt).
 * A threshold at/above v_N lets a rolling ball bank N slow frames and
 * freeze mid-slope, so the sleep threshold must be strictly below this. */
static float min_rolling_speed(float a) {
    float q = expf(-LIN_DAMP_PER_S * DT_S);
    float v = 0.0f;
    for (int k = 0; k < CONFIRM_FRAMES; k++)
        v = q * (v + a * DT_S);
    return v;
}

/* Lower bound: with nonzero restitution a resting contact re-emits a normal
 * speed of about e*g*dt every frame (then damping q). Under this a body
 * cannot bank N slow frames, so settled balls never sleep. */
#define E_GROUND 0.10f
static float jitter_floor(void) {
    float q = expf(-LIN_DAMP_PER_S * DT_S);
    return q * E_GROUND * G_MS2 * DT_S;
}

/* ============================== scene model ============================ */

#define BALL_R        0.14f
#define SENT_N        8
#define BED_N         30
#define SLOPE_X       12.0f   /* horizontal length of slope, m */
#define SENT_WORLD_XMAX (SLOPE_X + 14.0f)
#define BED_WORLD_XMAX  12.0f
#define RUN_FRAMES    1800    /* 30 s, lets the shallow-grade crawl settle */
#define WAKE_SPEED    0.10f   /* contact closing speed that wakes a sleeper */
#define E_BALL        0.20f
#define REST_SPEED    0.40f   /* slower closing contacts resolve inelastic */
#define WALL_E        0.30f
#define MU_BALL       0.25f
#define ORACLE_W      120     /* 2 s look-ahead for the no-sleep oracle */
#define ORACLE_DELTA  0.25f   /* m of oracle motion labelling a bad freeze */

typedef struct { float deg, sn, cs, tn; } Terrain;

static Terrain terrain_make(float deg) {
    Terrain t;
    t.deg = deg;
    float th = deg * PI_F / 180.0f;
    t.sn = sinf(th); t.cs = cosf(th); t.tn = tanf(th);
    return t;
}

typedef struct {
    Terrain t;
    float slope_x;         /* x < slope_x is sloped ground */
    float wall_max;        /* right wall plane (left wall at x=0) */
    float open_grade_xmax; /* sleep x below this is a mid-slope false sleep */
    int collide;           /* run sphere-sphere contacts (off for probes) */
} Scene;

static float ground_height(const Scene* s, float x) {
    if (x < s->slope_x) return (s->slope_x - x) * s->t.tn;
    return 0.0f;
}

static V3 ground_normal(const Scene* s, float x) {
    V3 n;
    if (x < s->slope_x) { n.x = s->t.sn; n.y = s->t.cs; }
    else                { n.x = 0.0f;    n.y = 1.0f; }
    n.z = 0.0f;
    return n;
}

static V3 ground_tangent(const Scene* s, float x) {
    V3 q;
    if (x < s->slope_x) { q.x = s->t.cs; q.y = -s->t.sn; }
    else                { q.x = 1.0f;    q.y = 0.0f; }
    q.z = 0.0f;
    return q;
}

typedef struct {
    V3 p, v;
    float r;
    SleepState sl;
} Ball;

typedef struct {
    int ball, frame;
    int on_slope;
    int wake_frame;
    float oracle_disp;
    int harmful;
} SleepEvent;

typedef struct { SleepEvent* ev; int count, cap; } EventLog;

static void events_push(EventLog* l, SleepEvent e) {
    if (l->count >= l->cap) {
        l->cap = l->cap ? l->cap * 2 : 64;
        l->ev = realloc(l->ev, (size_t)l->cap * sizeof(*l->ev));
    }
    e.wake_frame = -1;
    e.oracle_disp = 0.0f;
    e.harmful = 0;
    l->ev[l->count++] = e;
}

static unsigned rng_state = 0x12345678u;
static float rand_signed(float amp) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 17;
    rng_state ^= rng_state << 5;
    float u = (rng_state >> 8) * (1.0f / 16777216.0f);
    return (u * 2.0f - 1.0f) * amp;
}

/* Widely spaced, zero-interaction rolling train on the grade. */
static Scene make_sentinel_scene(float deg) {
    Scene s;
    s.t = terrain_make(deg);
    s.slope_x = SLOPE_X;
    s.wall_max = SENT_WORLD_XMAX;
    s.open_grade_xmax = SLOPE_X - 1.0f;
    s.collide = 0;
    return s;
}

static void build_sentinels(Ball* b, int n, float deg) {
    Scene s = make_sentinel_scene(deg);
    for (int i = 0; i < n; i++) {
        Ball* q = &b[i];
        memset(q, 0, sizeof(*q));
        q->r = BALL_R;
        float x = (n == 1) ? 1.0f : 0.60f + (float)i * 1.20f;
        q->p.x = x;
        q->p.y = ground_height(&s, x) + BALL_R + 0.002f;
        sleep_init(&q->sl);
    }
}

/* Compact agitated pile on a flat enclosed floor. */
static Scene make_bed_scene(void) {
    Scene s;
    s.t = terrain_make(0.0f);
    s.slope_x = 0.0f;
    s.wall_max = BED_WORLD_XMAX;
    s.open_grade_xmax = -1.0f;
    s.collide = 1;
    return s;
}

/* Agitated monolayer on a flat floor: lots of horizontal contacts, no tall
 * stack for a single-pass solver to churn on. */
static void build_bed(Ball* b, int n, unsigned seed) {
    rng_state = seed ? seed : 0x5eed1234u;
    const float sp = 2.0f * BALL_R + 0.08f; /* 0.36 */
    const float x0 = 0.55f;
    for (int i = 0; i < n; i++) {
        Ball* q = &b[i];
        memset(q, 0, sizeof(*q));
        q->r = BALL_R;
        q->p.x = x0 + (float)i * sp;
        q->p.y = BALL_R + 0.01f + ((i & 1) ? 0.02f : 0.0f);
        q->v.x = rand_signed(1.60f);
        q->v.y = rand_signed(0.20f);
        sleep_init(&q->sl);
    }
}

typedef struct {
    long pair_total, pair_skip, asleep_frames;
    int  sleep_events, false_slope;
} StepCounters;

static void wake_ball(Ball* b, int idx, int frame, EventLog* log) {
    if (!b[idx].sl.asleep) return;
    sleep_wake(&b[idx].sl);
    if (log) {
        for (int k = log->count - 1; k >= 0; k--) {
            if (log->ev[k].ball == idx && log->ev[k].wake_frame < 0) {
                log->ev[k].wake_frame = frame;
                break;
            }
        }
    }
}

static void world_step(Ball* b, int n, const Scene* sc,
                       const SleepCfg* cfg, int use_sleep,
                       int frame, EventLog* log, StepCounters* cnt) {
    const float damp = expf(-LIN_DAMP_PER_S * DT_S);

    for (int i = 0; i < n; i++) {
        if (use_sleep && b[i].sl.asleep) continue;

        b[i].v.y -= G_MS2 * DT_S;
        b[i].p = vadd(b[i].p, vscale(b[i].v, DT_S));

        if (b[i].p.x < b[i].r) {
            b[i].p.x = b[i].r;
            if (b[i].v.x < 0) b[i].v.x = -WALL_E * b[i].v.x;
        }
        if (b[i].p.x > sc->wall_max - b[i].r) {
            b[i].p.x = sc->wall_max - b[i].r;
            if (b[i].v.x > 0) b[i].v.x = -WALL_E * b[i].v.x;
        }

        float gh = ground_height(sc, b[i].p.x);
        float pen = (gh + b[i].r) - b[i].p.y;
        if (pen > 0.0f) {
            b[i].p.y += pen;
            V3 nrm = ground_normal(sc, b[i].p.x);
            V3 td  = ground_tangent(sc, b[i].p.x);

            float vn = vdot(b[i].v, nrm);
            if (vn < 0.0f)
                b[i].v = vsub(b[i].v, vscale(nrm, (1.0f + E_GROUND) * vn));

            int on_grade = b[i].p.x < sc->slope_x;
            float sin_here = on_grade ? sc->t.sn : 0.0f;
            float cos_here = on_grade ? sc->t.cs : 1.0f;

            /* particle gravity kick -> solid-sphere rolling (keep 5/7) */
            float vt = vdot(b[i].v, td);
            vt -= (1.0f - SOLID_SPHERE_K) * G_MS2 * sin_here * DT_S;

            float resist = ROLL_RESIST * G_MS2 * cos_here * DT_S;
            if (vt > 0.0f)      vt = fmaxf(0.0f, vt - resist);
            else if (vt < 0.0f) vt = fminf(0.0f, vt + resist);

            float vn2 = vdot(b[i].v, nrm);
            b[i].v = vadd(vscale(nrm, vn2), vscale(td, vt));
        }

        b[i].v = vscale(b[i].v, damp);
    }

    /* sphere-sphere contacts (disabled for isolated rolling probes) */
    if (sc->collide) for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            cnt->pair_total++;
            if (use_sleep && b[i].sl.asleep && b[j].sl.asleep) {
                cnt->pair_skip++;
                continue;
            }

            float dx = b[j].p.x - b[i].p.x;
            float dy = b[j].p.y - b[i].p.y;
            float d2 = dx*dx + dy*dy;
            float rr = b[i].r + b[j].r;
            if (d2 >= rr * rr) continue;

            float d = sqrtf(d2);
            if (d < 1e-8f) d = 1e-8f;
            float nx = dx / d, ny = dy / d;
            V3 nrm = { nx, ny, 0.0f };

            V3 rv = vsub(b[j].v, b[i].v);
            float vn = vdot(rv, nrm);

            if (use_sleep && vn < -WAKE_SPEED) {
                wake_ball(b, i, frame, log);
                wake_ball(b, j, frame, log);
            }

            float imi = (use_sleep && b[i].sl.asleep) ? 0.0f : 1.0f;
            float imj = (use_sleep && b[j].sl.asleep) ? 0.0f : 1.0f;
            float imsum = imi + imj;
            if (imsum <= 0.0f) continue;

            float overlap = rr - d;
            float pc = overlap / imsum;
            b[i].p.x -= nx * pc * imi;
            b[i].p.y -= ny * pc * imi;
            b[j].p.x += nx * pc * imj;
            b[j].p.y += ny * pc * imj;

            /* restitution only for genuine impacts; resting chatter is e=0 */
            float eeff = (-vn > REST_SPEED) ? E_BALL : 0.0f;
            float jimp = -(1.0f + eeff) * vn / imsum;
            if (jimp > 0.0f) {
                b[i].v = vsub(b[i].v, vscale(nrm, jimp * imi));
                b[j].v = vadd(b[j].v, vscale(nrm, jimp * imj));

                rv = vsub(b[j].v, b[i].v);
                V3 tng = { -ny, nx, 0.0f };
                float jt = -vdot(rv, tng) / imsum;
                float maxjt = MU_BALL * jimp;
                if (jt >  maxjt) jt =  maxjt;
                if (jt < -maxjt) jt = -maxjt;
                b[i].v = vsub(b[i].v, vscale(tng, jt * imi));
                b[j].v = vadd(b[j].v, vscale(tng, jt * imj));
            }
        }
    }

    if (use_sleep) {
        for (int i = 0; i < n; i++) {
            if (sleep_update(&b[i].sl, &b[i].v, cfg)) {
                cnt->sleep_events++;
                int on_open_grade = b[i].p.x < sc->open_grade_xmax;
                if (on_open_grade) cnt->false_slope++;
                if (log) {
                    SleepEvent e;
                    memset(&e, 0, sizeof(e));
                    e.ball = i;
                    e.frame = frame;
                    e.on_slope = on_open_grade;
                    events_push(log, e);
                }
            }
            if (b[i].sl.asleep) cnt->asleep_frames++;
        }
    }
}

/* ============================== run results =========================== */

typedef struct {
    int    n;
    int    false_slope;
    int    false_oracle;
    int    sleep_events;
    int    asleep_end;
    long   asleep_frames;
    long   pair_total, pair_skip;
    double ratio_end;
    double duty_ratio;
    double pair_skip_ratio;
} Result;

static float* run_oracle(const Scene* sc, int n, unsigned seed, int is_bed) {
    const int total = RUN_FRAMES + ORACLE_W;
    Ball* b = calloc((size_t)n, sizeof(Ball));
    float* trace = malloc((size_t)(total + 1) * n * 2 * sizeof(float));
    if (is_bed) build_bed(b, n, seed); else build_sentinels(b, n, sc->t.deg);

    SleepCfg off = { 1e30f, 1e30f, CONFIRM_FRAMES };
    for (int i = 0; i < n; i++) {
        trace[i * 2 + 0] = b[i].p.x;
        trace[i * 2 + 1] = b[i].p.y;
    }
    StepCounters z; memset(&z, 0, sizeof(z));
    for (int f = 1; f <= total; f++) {
        world_step(b, n, sc, &off, 0, f, NULL, &z);
        for (int i = 0; i < n; i++) {
            trace[(f * n + i) * 2 + 0] = b[i].p.x;
            trace[(f * n + i) * 2 + 1] = b[i].p.y;
        }
    }
    free(b);
    return trace;
}

static Result run_scene(const Scene* sc, int n, float threshold,
                        unsigned seed, int is_bed,
                        int use_oracle, int verbose) {
    Ball* b = calloc((size_t)n, sizeof(Ball));
    if (is_bed) build_bed(b, n, seed); else build_sentinels(b, n, sc->t.deg);

    SleepCfg cfg = { threshold, threshold, CONFIRM_FRAMES };
    EventLog log = { 0 };
    StepCounters cnt; memset(&cnt, 0, sizeof(cnt));

    for (int f = 1; f <= RUN_FRAMES; f++)
        world_step(b, n, sc, &cfg, 1, f, &log, &cnt);

    Result r; memset(&r, 0, sizeof(r));
    r.n = n;
    r.false_slope = cnt.false_slope;
    r.sleep_events = cnt.sleep_events;
    r.asleep_frames = cnt.asleep_frames;
    r.pair_total = cnt.pair_total;
    r.pair_skip = cnt.pair_skip;
    for (int i = 0; i < n; i++) if (b[i].sl.asleep) r.asleep_end++;
    r.ratio_end = (double)r.asleep_end / n;
    r.duty_ratio = (double)cnt.asleep_frames / ((double)n * RUN_FRAMES);
    r.pair_skip_ratio = cnt.pair_total ? (double)cnt.pair_skip / cnt.pair_total : 0.0;

    if (use_oracle) {
        float* tr = run_oracle(sc, n, seed, is_bed);
        int bad = 0;
        for (int k = 0; k < log.count; k++) {
            SleepEvent* e = &log.ev[k];
            int f2 = e->frame + ORACLE_W;
            int i = e->ball;
            float x0 = tr[(e->frame * n + i) * 2 + 0];
            float y0 = tr[(e->frame * n + i) * 2 + 1];
            float x1 = tr[(f2 * n + i) * 2 + 0];
            float y1 = tr[(f2 * n + i) * 2 + 1];
            e->oracle_disp = sqrtf((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
            int stayed = (e->wake_frame < 0) ||
                         (e->wake_frame - e->frame > ORACLE_W);
            e->harmful = (e->oracle_disp > ORACLE_DELTA) && stayed;
            if (e->harmful) bad++;
            if (verbose)
                fprintf(stderr, "  ev ball=%d f=%d open_grade=%d wake=%d "
                                "disp=%.3f harmful=%d\n",
                        e->ball, e->frame, e->on_slope, e->wake_frame,
                        e->oracle_disp, e->harmful);
        }
        r.false_oracle = bad;
        free(tr);
    }

    free(log.ev);
    free(b);
    return r;
}

/* Smallest threshold that puts an isolated from-rest ball to sleep on the
 * open grade, found by bisection: the simulated false-sleep knee. */
static float simulated_knee(float deg) {
    float lo = 0.001f, hi = 3.0f;
    Scene sc = make_sentinel_scene(deg);
    for (int it = 0; it < 26; it++) {
        float mid = 0.5f * (lo + hi);
        Ball b;
        build_sentinels(&b, 1, deg);
        SleepCfg cfg = { mid, mid, CONFIRM_FRAMES };
        StepCounters cnt; memset(&cnt, 0, sizeof(cnt));
        for (int f = 1; f <= 240; f++) {
            world_step(&b, 1, &sc, &cfg, 1, f, NULL, &cnt);
            if (cnt.false_slope > 0) break;
        }
        if (cnt.false_slope > 0) hi = mid; else lo = mid;
    }
    return hi;
}

/* ============================== self tests ============================ */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static int run_self_tests(void) {
    SleepCfg cfg = { 0.05f, 0.05f, 30 };

    V3 still = { 0.01f, -0.01f, 0.0f };
    check("slow_both_components", sleep_is_slow(&still, &cfg) == 1);

    /* THE anti-false-sleep case: zero vertical speed, horizontal motion. */
    V3 horizontal = { 2.0f, 0.0f, 0.0f };
    check("vy_zero_but_horizontal_not_slow", sleep_is_slow(&horizontal, &cfg) == 0);
    SleepState s1; sleep_init(&s1);
    int slept_h = 0;
    for (int f = 0; f < 120; f++)
        slept_h |= sleep_update(&s1, &horizontal, &cfg);
    check("horizontal_mover_never_sleeps", slept_h == 0 && s1.asleep == 0);

    V3 vertical = { 0.0f, -1.5f, 0.0f };
    check("vertical_motion_not_slow", sleep_is_slow(&vertical, &cfg) == 0);

    V3 zdepth = { 0.0f, 0.0f, 0.4f };
    check("z_horizontal_component_gated", sleep_is_slow(&zdepth, &cfg) == 0);

    SleepState s2; sleep_init(&s2);
    V3 tiny = { 0.01f, 0.01f, 0.0f };
    int transition = 0;
    for (int f = 1; f <= CONFIRM_FRAMES; f++) {
        int e = sleep_update(&s2, &tiny, &cfg);
        if (f < CONFIRM_FRAMES) { if (e) transition = -1; }
        else transition = e;
    }
    check("confirm_frames_then_sleep", transition == 1 && s2.asleep == 1);

    SleepState s3; sleep_init(&s3);
    for (int f = 0; f < CONFIRM_FRAMES - 1; f++) sleep_update(&s3, &tiny, &cfg);
    V3 fast = { 0.3f, 0.0f, 0.0f };
    sleep_update(&s3, &fast, &cfg);
    check("fast_frame_resets_counter", s3.slow_count == 0 && s3.asleep == 0);

    sleep_wake(&s2);
    check("wake_clears_state", s2.asleep == 0 && s2.slow_count == 0);

    V3 edge = { 0.05f, 0.0f, 0.0f };
    check("threshold_boundary_strict", sleep_is_slow(&edge, &cfg) == 0);

    float a_min = slope_accel(GRADES_DEG[0]);
    float v_min = min_rolling_speed(a_min);
    check("shallow_grade_still_rolls", a_min > 0.05f);
    check("threshold_below_rolling_min", THRESHOLD_SELECTED < v_min);
    check("threshold_above_jitter_floor", THRESHOLD_SELECTED > jitter_floor());

    SleepCfg sel = { THRESHOLD_SELECTED, THRESHOLD_SELECTED, CONFIRM_FRAMES };
    V3 hcase = { THRESHOLD_SELECTED + 0.001f, 0.0f, 0.0f };
    check("selected_gate_blocks_horizontal", sleep_is_slow(&hcase, &sel) == 0);

    /* acceptance: isolated sentinels on all three grades have zero false
     * sleeps (geometry + oracle) at the selected threshold. */
    int all_zero = 1;
    for (int g = 0; g < GRADE_COUNT; g++) {
        Scene sc = make_sentinel_scene(GRADES_DEG[g]);
        Result r = run_scene(&sc, SENT_N, THRESHOLD_SELECTED, 1u, 0, 1, 0);
        if (r.false_slope != 0 || r.false_oracle != 0) all_zero = 0;
    }
    check("grades_zero_false_sleep", all_zero);

    /* settling pile must put a clear majority to sleep at the threshold */
    Scene bed = make_bed_scene();
    Result rb = run_scene(&bed, BED_N, THRESHOLD_SELECTED, 0xBED, 1, 0, 0);
    check("bed_majority_asleep", rb.ratio_end >= 0.8);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}

/* ============================== CLI modes ============================== */

static int cmd_derive(void) {
    printf("g=%.2f dt=%.6f confirm_frames=%d window=%.3fs\n",
           G_MS2, DT_S, CONFIRM_FRAMES, CONFIRM_FRAMES * DT_S);
    printf("solid_sphere_factor=%.5f rolling_resistance=%.3f damping=%.3f\n",
           SOLID_SPHERE_K, ROLL_RESIST, LIN_DAMP_PER_S);
    printf("jitter_floor=%.4f\n", jitter_floor());
    for (int g = 0; g < GRADE_COUNT; g++) {
        float a = slope_accel(GRADES_DEG[g]);
        printf("grade_deg=%.0f a_net=%.4f discrete_vN=%.4f simulated_knee=%.4f\n",
               GRADES_DEG[g], a, min_rolling_speed(a),
               simulated_knee(GRADES_DEG[g]));
    }
    printf("threshold_selected=%.4f\n", THRESHOLD_SELECTED);
    return 0;
}

static int cmd_knee_csv(void) {
    printf("grade_deg,a_net,v_min_window,simulated_knee\n");
    for (int g = 0; g < GRADE_COUNT; g++) {
        float a = slope_accel(GRADES_DEG[g]);
        printf("%.0f,%.5f,%.5f,%.5f\n", GRADES_DEG[g], a,
               min_rolling_speed(a), simulated_knee(GRADES_DEG[g]));
    }
    printf("jitter_floor,%.5f\n", jitter_floor());
    return 0;
}

static int cmd_sweep(const char* path) {
    FILE* f = fopen(path, "w");
    if (!f) { perror("fopen sweep"); return 2; }
    fprintf(f, "threshold,false_5,false_10,false_20,"
               "bed_asleep,bed_n,bed_ratio,bed_duty,bed_skip,"
               "total_ratio_5,total_ratio_10,total_ratio_20\n");

    const int NSTEPS = 46;
    const unsigned BEDSEED = 0xBEDu;
    for (int k = 0; k < NSTEPS; k++) {
        float u = (float)k / (NSTEPS - 1);
        float thr = 0.008f * powf(1.6f / 0.008f, u); /* log grid */

        Scene bedsc = make_bed_scene();
        Result bed = run_scene(&bedsc, BED_N, thr, BEDSEED, 1, 0, 0);

        int fs[3]; double tr[3];
        for (int g = 0; g < GRADE_COUNT; g++) {
            Scene sc = make_sentinel_scene(GRADES_DEG[g]);
            Result s = run_scene(&sc, SENT_N, thr, 1u, 0, 0, 0);
            fs[g] = s.false_slope;
            tr[g] = (double)(s.asleep_end + bed.asleep_end) / (SENT_N + BED_N);
        }
        fprintf(f, "%.5f,%d,%d,%d,%d,%d,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n",
                thr, fs[0], fs[1], fs[2],
                bed.asleep_end, BED_N, bed.ratio_end, bed.duty_ratio,
                bed.pair_skip_ratio, tr[0], tr[1], tr[2]);
        printf("threshold=%.4f false=(%d,%d,%d) bed_ratio=%.2f\n",
               thr, fs[0], fs[1], fs[2], bed.ratio_end);
    }
    fclose(f);
    printf("wrote %s\n", path);
    return 0;
}

static int cmd_grades_json(float thr) {
    printf("{\n");
    printf("  \"threshold\": %.4f,\n", thr);
    printf("  \"confirm_frames\": %d,\n", CONFIRM_FRAMES);
    printf("  \"dt\": %.6f,\n", DT_S);
    printf("  \"run_seconds\": %.1f,\n", RUN_FRAMES * DT_S);
    printf("  \"sentinel_n\": %d, \"bed_n\": %d,\n", SENT_N, BED_N);
    printf("  \"jitter_floor\": %.4f,\n", jitter_floor());
    printf("  \"grades\": [\n");
    for (int g = 0; g < GRADE_COUNT; g++) {
        float deg = GRADES_DEG[g];
        float a = slope_accel(deg), vm = min_rolling_speed(a);
        /* one self-contained run per grade: slope sentinels + settling pile */
        Scene sc = make_sentinel_scene(deg);
        Result s = run_scene(&sc, SENT_N, thr, 1u, 0, 1, 1);
        Scene bedsc = make_bed_scene();
        Result bed = run_scene(&bedsc, BED_N, thr, 0xBEDu + (unsigned)g, 1, 0, 0);

        long pt = s.pair_total + bed.pair_total;
        long ps = s.pair_skip + bed.pair_skip;
        double total_ratio =
            (double)(s.asleep_end + bed.asleep_end) / (SENT_N + BED_N);
        double duty =
            (double)(s.asleep_frames + bed.asleep_frames) /
            ((double)(SENT_N + BED_N) * RUN_FRAMES);

        printf("    {\"grade_deg\": %.0f, \"a_net\": %.4f, "
               "\"v_min_window\": %.4f, \"false_sleeps\": %d, "
               "\"false_oracle\": %d, \"sentinel_asleep\": %d, "
               "\"bed_sleep_ratio\": %.4f, \"sleep_ratio\": %.4f, "
               "\"duty_ratio\": %.4f, \"pair_skip_ratio\": %.4f}%s\n",
               deg, a, vm, s.false_slope, s.false_oracle, s.asleep_end,
               bed.ratio_end, total_ratio, duty,
               pt ? (double)ps / pt : 0.0,
               g + 1 < GRADE_COUNT ? "," : "");
    }
    printf("  ]\n}\n");
    return 0;
}

static int cmp_double(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return (x > y) - (x < y);
}

static double bench_one(const Scene* sc, int n, int frames, float thr,
                        int use_sleep, unsigned seed,
                        long* skip, long* total, int* nasleep) {
    Ball* b = calloc((size_t)n, sizeof(Ball));
    build_bed(b, n, seed);
    SleepCfg cfg = { thr, thr, CONFIRM_FRAMES };
    StepCounters cnt; memset(&cnt, 0, sizeof(cnt));
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (int f = 1; f <= frames; f++)
        world_step(b, n, sc, &cfg, use_sleep, f, NULL, &cnt);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double ms = (t1.tv_sec - t0.tv_sec) * 1e3 +
                (t1.tv_nsec - t0.tv_nsec) / 1e6;
    if (skip) *skip = cnt.pair_skip;
    if (total) *total = cnt.pair_total;
    if (nasleep) {
        *nasleep = 0;
        for (int i = 0; i < n; i++) *nasleep += b[i].sl.asleep;
    }
    free(b);
    return ms;
}

static int cmd_bench(void) {
    const int n = 200, frames = 2400, reps = 9;
    Scene sc = make_bed_scene();
    sc.wall_max = 0.55f + (float)n * (2.0f * BALL_R + 0.08f) + 1.0f;

    double on[16], off[16];
    long skip_last = 0, total_last = 1;
    int asleep_last = 0;
    for (int rep = 0; rep < reps; rep++) {
        unsigned seed = 777u + (unsigned)rep;
        /* interleave and alternate order to spread scheduling/thermal noise */
        if (rep & 1) {
            off[rep] = bench_one(&sc, n, frames, 1e30f, 0, seed, 0, 0, 0);
            on[rep]  = bench_one(&sc, n, frames, THRESHOLD_SELECTED, 1, seed,
                                 &skip_last, &total_last, &asleep_last);
        } else {
            on[rep]  = bench_one(&sc, n, frames, THRESHOLD_SELECTED, 1, seed,
                                 &skip_last, &total_last, &asleep_last);
            off[rep] = bench_one(&sc, n, frames, 1e30f, 0, seed, 0, 0, 0);
        }
    }
    qsort(on, reps, sizeof(double), cmp_double);
    qsort(off, reps, sizeof(double), cmp_double);
    double med_on = on[reps / 2], med_off = off[reps / 2];
    printf("balls=%d frames=%d reps=%d sleep_on_median_ms=%.2f "
           "sleep_off_median_ms=%.2f median_speedup=%.2fx "
           "best_speedup=%.2fx end_asleep=%d pair_skip_ratio=%.3f\n",
           n, frames, reps, med_on, med_off, med_off / med_on,
           off[0] / on[0], asleep_last,
           (double)skip_last / total_last);
    return 0;
}

static int cmd_trace(const char* which, float thr) {
    if (!strcmp(which, "bed")) {
        Scene sc = make_bed_scene();
        Ball* b = calloc(BED_N, sizeof(Ball));
        build_bed(b, BED_N, 0xBEDu);
        SleepCfg cfg = { thr, thr, CONFIRM_FRAMES };
        StepCounters cnt; memset(&cnt, 0, sizeof(cnt));
        for (int f = 1; f <= RUN_FRAMES; f++) {
            world_step(b, BED_N, &sc, &cfg, 1, f, NULL, &cnt);
            if (f % 60 == 0 || f == 1) {
                int asleep = 0; float maxv = 0;
                for (int i = 0; i < BED_N; i++) {
                    if (b[i].sl.asleep) asleep++;
                    float sp = sqrtf(vdot(b[i].v, b[i].v));
                    if (sp > maxv) maxv = sp;
                }
                printf("f=%4d t=%5.2f asleep=%2d maxv=%.3f events=%d\n",
                       f, f * DT_S, asleep, maxv, cnt.sleep_events);
            }
        }
        for (int i = 0; i < BED_N; i++) {
            if (!b[i].sl.asleep)
                printf("  awake i=%2d p=(%.2f,%.2f) v=(%.3f,%.3f) sp=%.3f\n",
                       i, b[i].p.x, b[i].p.y, b[i].v.x, b[i].v.y,
                       sqrtf(vdot(b[i].v, b[i].v)));
        }
        free(b);
    } else {
        float deg = (float)atof(which);
        Scene sc = make_sentinel_scene(deg);
        Ball b[SENT_N];
        build_sentinels(b, SENT_N, deg);
        SleepCfg cfg = { thr, thr, CONFIRM_FRAMES };
        StepCounters cnt; memset(&cnt, 0, sizeof(cnt));
        for (int f = 1; f <= RUN_FRAMES; f++) {
            world_step(b, SENT_N, &sc, &cfg, 1, f, NULL, &cnt);
            if (f % 120 == 0) {
                int asleep = 0;
                for (int i = 0; i < SENT_N; i++) asleep += b[i].sl.asleep;
                printf("f=%4d asleep=%d false=%d\n", f, asleep, cnt.false_slope);
            }
        }
        for (int i = 0; i < SENT_N; i++)
            if (!b[i].sl.asleep)
                printf("  awake sentinel i=%d p=(%.2f,%.2f) v=(%.4f,%.4f) sp=%.4f count=%d\n",
                       i, b[i].p.x, b[i].p.y, b[i].v.x, b[i].v.y,
                       sqrtf(vdot(b[i].v, b[i].v)), b[i].sl.slow_count);
    }
    return 0;
}

int main(int argc, char** argv) {
    if (argc == 1) return run_self_tests();
    if (!strcmp(argv[1], "test"))   return run_self_tests();
    if (!strcmp(argv[1], "derive")) return cmd_derive();
    if (!strcmp(argv[1], "knee"))   return cmd_knee_csv();
    if (!strcmp(argv[1], "bench"))  return cmd_bench();
    if (!strcmp(argv[1], "sweep") && argc >= 3) return cmd_sweep(argv[2]);
    if (!strcmp(argv[1], "grades")) {
        float thr = (argc >= 3) ? (float)atof(argv[2]) : THRESHOLD_SELECTED;
        return cmd_grades_json(thr);
    }
    if (!strcmp(argv[1], "trace") && argc >= 4)
        return cmd_trace(argv[2], (float)atof(argv[3]));
    fprintf(stderr,
        "usage: %s [test|derive|knee|bench|sweep out.csv|grades T|trace bed|<deg> T]\n",
        argv[0]);
    return 2;
}
