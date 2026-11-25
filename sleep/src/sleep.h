#ifndef MINI_PHYS_SLEEP_H
#define MINI_PHYS_SLEEP_H

/*
 * Four-state ball sleep FSM.
 *
 * AWAKE    -- normal simulation.
 * DROWSY   -- first calm frames; sleep is provisional and any stimulus
 *             returns immediately to AWAKE.
 * ASLEEP   -- simulation is skipped until a collision or above-threshold
 *             activity bit requests wake-up.
 * WAKING   -- latched wake ramp; it cannot fall back towards ASLEEP before
 *             reaching AWAKE, which gives a finite wake-latency bound.
 */

#ifdef __cplusplus
extern "C" {
#endif

#define SLEEP_FRAME_HZ 60u

/* Calibration chosen by artifacts/sleep_tradeoff.png. */
#define SLEEP_SELECTED_DROWSY_FRAMES 6u
#define SLEEP_SELECTED_WAKE_FRAMES   2u

/* Replay calibration target and latency budget. */
#define SLEEP_FALSE_SLEEP_BUDGET 10u
#define SLEEP_WAKE_BUDGET_FRAMES 3u

typedef enum SleepState {
    SLEEP_AWAKE = 0,
    SLEEP_DROWSY = 1,
    SLEEP_ASLEEP = 2,
    SLEEP_WAKING = 3
} SleepState;

typedef enum SleepInputKind {
    SLEEP_INPUT_CALM = 0,
    SLEEP_INPUT_ACTIVE = 1,
    SLEEP_INPUT_HIT = 2
} SleepInputKind;

typedef struct SleepInput {
    unsigned int hit;     /* collision/impact this frame */
    unsigned int active;  /* speed/energy above the awake threshold */
} SleepInput;

typedef enum SleepEdgeId {
    SLEEP_EDGE_NONE = 0,
    SLEEP_EDGE_AWAKE_HOLD,
    SLEEP_EDGE_AWAKE_TO_DROWSY,
    SLEEP_EDGE_DROWSY_HOLD,
    SLEEP_EDGE_DROWSY_TO_AWAKE,
    SLEEP_EDGE_DROWSY_TO_ASLEEP,
    SLEEP_EDGE_ASLEEP_HOLD,
    SLEEP_EDGE_ASLEEP_TO_WAKING,
    SLEEP_EDGE_WAKING_HOLD,
    SLEEP_EDGE_WAKING_TO_AWAKE,
    SLEEP_EDGE_RECOVER_INVALID,
    SLEEP_EDGE_COUNT
} SleepEdgeId;

typedef struct SleepConfig {
    unsigned int drowsy_frames; /* D: calm frames required inside DROWSY */
    unsigned int wake_frames;   /* W: frames required inside WAKING */
} SleepConfig;

typedef struct SleepFsm {
    SleepState state;
    SleepConfig config;
    unsigned int drowsy_left;
    unsigned int wake_left;
    unsigned long long frame;
    SleepEdgeId last_edge;
    SleepInputKind last_input;
} SleepFsm;

int sleep_config_valid(SleepConfig config);
SleepConfig sleep_selected_config(void);
void sleep_init(SleepFsm* fsm, SleepConfig config);

/* hit has priority over active; both absent means calm. */
SleepInputKind sleep_input_kind(SleepInput input);

/* Returns 0 on a normal step, -1 for a null/invalid argument.  An invalid
 * encoded state is defensively recovered to AWAKE and reported through
 * last_edge, so corruption cannot persist as an unnamed fifth state. */
int sleep_step(SleepFsm* fsm, SleepInput input);

const char* sleep_state_name(SleepState state);
const char* sleep_edge_name(SleepEdgeId edge);

/* Graph edge budget in step() calls/frames.  D/W edges are config-sized. */
unsigned int sleep_edge_frame_budget(SleepEdgeId edge, SleepConfig config);

/* Edge-length and path bounds, in step() calls/frames. */
unsigned int sleep_awake_to_drowsy_bound(void);
unsigned int sleep_drowsy_to_asleep_bound(SleepConfig config);
unsigned int sleep_first_calm_to_asleep_bound(SleepConfig config);
unsigned int sleep_asleep_to_waking_bound(void);
unsigned int sleep_waking_to_awake_bound(SleepConfig config);
unsigned int sleep_hit_to_awake_bound(SleepConfig config);

/*
 * Wake-path monitor.  Begin immediately before submitting the collision
 * frame, then call observe() after every subsequent sleep_step().  It
 * remembers the first edge that fails to finish within its graph budget.
 * The first observed AWAKE freezes actual_frames; later observe() calls are
 * ignored, so it is safe to keep feeding the monitor from the main loop.
 */
typedef struct SleepWakeReport {
    int finished;
    int within_bound;
    unsigned int actual_frames;
    unsigned int bound_frames;
    SleepEdgeId stuck_edge;
    SleepState terminal_state;
} SleepWakeReport;

typedef struct SleepWakeWatch {
    SleepConfig config;
    unsigned int frames;
    unsigned int waking_frames;
    unsigned int completion_frame;
    SleepEdgeId stuck_edge;
    int violation_latched;
    int completed; /* latched on the first observed AWAKE */
} SleepWakeWatch;

void sleep_wake_watch_begin(SleepWakeWatch* watch, SleepConfig config);
void sleep_wake_watch_observe(SleepWakeWatch* watch, SleepState observed);
SleepWakeReport sleep_wake_watch_report(const SleepWakeWatch* watch,
                                        SleepState terminal_state);

/* Run the canonical ASLEEP --hit--> WAKING --> AWAKE scenario. */
SleepWakeReport sleep_run_canonical_wake(SleepConfig config);

#ifdef __cplusplus
}
#endif

#endif
