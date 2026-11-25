/*
 * sleep: four-state ball sleep FSM, completeness/transition self-tests,
 * wake-latency diagnostics, and the threshold sweep used for calibration.
 *
 * `make test`  runs the self-checks.
 * `make sweep` writes artifacts/sleep_tradeoff.csv.
 */
#include "sleep.h"

#include <stdio.h>
#include <string.h>

#define SWEEP_MAX_N 16u
#define FAULT_DEMONSTRATION_EXTRA 5u

static int failures = 0;

static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int sleep_config_valid(SleepConfig config) {
    return config.drowsy_frames >= 1u && config.drowsy_frames <= 100000u &&
           config.wake_frames >= 1u && config.wake_frames <= 100000u;
}

SleepConfig sleep_selected_config(void) {
    SleepConfig config = {
        SLEEP_SELECTED_DROWSY_FRAMES,
        SLEEP_SELECTED_WAKE_FRAMES
    };
    return config;
}

void sleep_init(SleepFsm* fsm, SleepConfig config) {
    if (fsm == 0) return;
    memset(fsm, 0, sizeof(*fsm));
    fsm->state = SLEEP_AWAKE;
    fsm->config = sleep_config_valid(config) ? config : sleep_selected_config();
    fsm->last_edge = SLEEP_EDGE_NONE;
    fsm->last_input = SLEEP_INPUT_CALM;
}

SleepInputKind sleep_input_kind(SleepInput input) {
    if (input.hit) return SLEEP_INPUT_HIT;
    if (input.active) return SLEEP_INPUT_ACTIVE;
    return SLEEP_INPUT_CALM;
}

int sleep_step(SleepFsm* fsm, SleepInput input) {
    if (fsm == 0) return -1;

    /* Defensive totality for a corrupted encoded state.  This runs before
     * config validation so even a simultaneously bad config cannot leave an
     * out-of-range state persisted. */
    if ((unsigned int)fsm->state > (unsigned int)SLEEP_WAKING) {
        fsm->last_input = sleep_input_kind(input);
        fsm->state = SLEEP_AWAKE;
        fsm->drowsy_left = 0;
        fsm->wake_left = 0;
        fsm->last_edge = SLEEP_EDGE_RECOVER_INVALID;
        fsm->frame++;
        return 0;
    }

    if (!sleep_config_valid(fsm->config)) return -1;

    SleepInputKind kind = sleep_input_kind(input);
    fsm->last_input = kind;

    switch (fsm->state) {
    case SLEEP_AWAKE:
        if (kind == SLEEP_INPUT_CALM) {
            fsm->state = SLEEP_DROWSY;
            fsm->drowsy_left = fsm->config.drowsy_frames;
            fsm->last_edge = SLEEP_EDGE_AWAKE_TO_DROWSY;
        } else {
            fsm->last_edge = SLEEP_EDGE_AWAKE_HOLD;
        }
        break;

    case SLEEP_DROWSY:
        if (kind != SLEEP_INPUT_CALM) {
            fsm->state = SLEEP_AWAKE;
            fsm->drowsy_left = 0;
            fsm->last_edge = SLEEP_EDGE_DROWSY_TO_AWAKE;
        } else {
            if (fsm->drowsy_left == 0)
                fsm->drowsy_left = fsm->config.drowsy_frames;
            fsm->drowsy_left--;
            if (fsm->drowsy_left == 0) {
                fsm->state = SLEEP_ASLEEP;
                fsm->wake_left = 0;
                fsm->last_edge = SLEEP_EDGE_DROWSY_TO_ASLEEP;
            } else {
                fsm->last_edge = SLEEP_EDGE_DROWSY_HOLD;
            }
        }
        break;

    case SLEEP_ASLEEP:
        if (kind != SLEEP_INPUT_CALM) {
            fsm->state = SLEEP_WAKING;
            fsm->wake_left = fsm->config.wake_frames;
            fsm->last_edge = SLEEP_EDGE_ASLEEP_TO_WAKING;
        } else {
            fsm->last_edge = SLEEP_EDGE_ASLEEP_HOLD;
        }
        break;

    case SLEEP_WAKING:
        if (fsm->wake_left == 0)
            fsm->wake_left = fsm->config.wake_frames;
        /* The wake ramp is latched: further calm input cannot move the
         * machine back towards DROWSY/ASLEEP and restart the path. */
        fsm->wake_left--;
        if (fsm->wake_left == 0) {
            fsm->state = SLEEP_AWAKE;
            fsm->last_edge = SLEEP_EDGE_WAKING_TO_AWAKE;
        } else {
            fsm->last_edge = SLEEP_EDGE_WAKING_HOLD;
        }
        break;
    }

    fsm->frame++;
    return 0;
}

const char* sleep_state_name(SleepState state) {
    switch (state) {
    case SLEEP_AWAKE:  return "AWAKE";
    case SLEEP_DROWSY: return "DROWSY";
    case SLEEP_ASLEEP: return "ASLEEP";
    case SLEEP_WAKING: return "WAKING";
    }
    return "INVALID_STATE";
}

const char* sleep_edge_name(SleepEdgeId edge) {
    switch (edge) {
    case SLEEP_EDGE_NONE: return "NONE";
    case SLEEP_EDGE_AWAKE_HOLD: return "AWAKE_HOLD";
    case SLEEP_EDGE_AWAKE_TO_DROWSY: return "AWAKE->DROWSY";
    case SLEEP_EDGE_DROWSY_HOLD: return "DROWSY_HOLD";
    case SLEEP_EDGE_DROWSY_TO_AWAKE: return "DROWSY->AWAKE";
    case SLEEP_EDGE_DROWSY_TO_ASLEEP: return "DROWSY->ASLEEP";
    case SLEEP_EDGE_ASLEEP_HOLD: return "ASLEEP_HOLD";
    case SLEEP_EDGE_ASLEEP_TO_WAKING: return "ASLEEP->WAKING";
    case SLEEP_EDGE_WAKING_HOLD: return "WAKING_HOLD";
    case SLEEP_EDGE_WAKING_TO_AWAKE: return "WAKING->AWAKE";
    case SLEEP_EDGE_RECOVER_INVALID: return "INVALID->AWAKE";
    case SLEEP_EDGE_COUNT: break;
    }
    return "UNKNOWN_EDGE";
}

unsigned int sleep_edge_frame_budget(SleepEdgeId edge, SleepConfig config) {
    switch (edge) {
    case SLEEP_EDGE_AWAKE_TO_DROWSY:
    case SLEEP_EDGE_DROWSY_TO_AWAKE:
    case SLEEP_EDGE_ASLEEP_TO_WAKING:
        return 1u;
    case SLEEP_EDGE_DROWSY_TO_ASLEEP:
        return sleep_config_valid(config) ? config.drowsy_frames : 0u;
    case SLEEP_EDGE_WAKING_TO_AWAKE:
        return sleep_config_valid(config) ? config.wake_frames : 0u;
    case SLEEP_EDGE_AWAKE_HOLD:
    case SLEEP_EDGE_DROWSY_HOLD:
    case SLEEP_EDGE_ASLEEP_HOLD:
    case SLEEP_EDGE_WAKING_HOLD:
        return 1u;
    default:
        return 0u;
    }
}

unsigned int sleep_awake_to_drowsy_bound(void) { return 1u; }

unsigned int sleep_drowsy_to_asleep_bound(SleepConfig config) {
    return sleep_edge_frame_budget(SLEEP_EDGE_DROWSY_TO_ASLEEP, config);
}

unsigned int sleep_first_calm_to_asleep_bound(SleepConfig config) {
    if (!sleep_config_valid(config)) return 0u;
    return sleep_awake_to_drowsy_bound() +
           sleep_drowsy_to_asleep_bound(config);
}

unsigned int sleep_asleep_to_waking_bound(void) {
    return sleep_edge_frame_budget(SLEEP_EDGE_ASLEEP_TO_WAKING,
                                   sleep_selected_config());
}

unsigned int sleep_waking_to_awake_bound(SleepConfig config) {
    return sleep_edge_frame_budget(SLEEP_EDGE_WAKING_TO_AWAKE, config);
}

unsigned int sleep_hit_to_awake_bound(SleepConfig config) {
    if (!sleep_config_valid(config)) return 0u;
    return sleep_edge_frame_budget(SLEEP_EDGE_ASLEEP_TO_WAKING, config) +
           sleep_edge_frame_budget(SLEEP_EDGE_WAKING_TO_AWAKE, config);
}

void sleep_wake_watch_begin(SleepWakeWatch* watch, SleepConfig config) {
    if (watch == 0) return;
    memset(watch, 0, sizeof(*watch));
    watch->config = sleep_config_valid(config) ? config : sleep_selected_config();
    watch->stuck_edge = SLEEP_EDGE_NONE;
}

void sleep_wake_watch_observe(SleepWakeWatch* watch, SleepState observed) {
    if (watch == 0 || watch->completed) return;

    watch->frames++;
    if (observed == SLEEP_WAKING) watch->waking_frames++;

    unsigned int bound = sleep_hit_to_awake_bound(watch->config);

    if (!watch->violation_latched) {
        if (watch->frames == 1u &&
            observed != SLEEP_WAKING && observed != SLEEP_AWAKE) {
            watch->stuck_edge = SLEEP_EDGE_ASLEEP_TO_WAKING;
            watch->violation_latched = 1;
        } else if (observed == SLEEP_WAKING &&
                   watch->waking_frames > watch->config.wake_frames) {
            watch->stuck_edge = SLEEP_EDGE_WAKING_TO_AWAKE;
            watch->violation_latched = 1;
        } else if (watch->frames > bound) {
            watch->stuck_edge = (observed == SLEEP_ASLEEP)
                ? SLEEP_EDGE_ASLEEP_TO_WAKING
                : SLEEP_EDGE_WAKING_TO_AWAKE;
            watch->violation_latched = 1;
        }
    }

    if (observed == SLEEP_AWAKE) {
        watch->completed = 1;
        watch->completion_frame = watch->frames;
    }
}

SleepWakeReport sleep_wake_watch_report(const SleepWakeWatch* watch,
                                        SleepState terminal_state) {
    SleepWakeReport report;
    memset(&report, 0, sizeof(report));
    if (watch == 0) {
        report.stuck_edge = SLEEP_EDGE_NONE;
        return report;
    }

    report.bound_frames = sleep_hit_to_awake_bound(watch->config);
    report.actual_frames = watch->completed
        ? watch->completion_frame
        : watch->frames;
    report.finished = watch->completed || (terminal_state == SLEEP_AWAKE);
    report.terminal_state = terminal_state;
    report.stuck_edge = watch->stuck_edge;

    if (report.finished && report.actual_frames <= report.bound_frames &&
        report.stuck_edge == SLEEP_EDGE_NONE) {
        report.within_bound = 1;
    } else {
        report.within_bound = 0;
        if (report.stuck_edge == SLEEP_EDGE_NONE &&
            report.actual_frames > report.bound_frames) {
            report.stuck_edge = (terminal_state == SLEEP_ASLEEP)
                ? SLEEP_EDGE_ASLEEP_TO_WAKING
                : SLEEP_EDGE_WAKING_TO_AWAKE;
        }
    }
    return report;
}

SleepWakeReport sleep_run_canonical_wake(SleepConfig config) {
    SleepFsm fsm;
    sleep_init(&fsm, config);
    config = fsm.config;

    const SleepInput calm = {0u, 0u};
    const SleepInput hit = {1u, 0u};

    /* AWAKE --one calm frame--> DROWSY --D calm frames--> ASLEEP. */
    (void)sleep_step(&fsm, calm);
    for (unsigned int i = 0;
         i < config.drowsy_frames && fsm.state != SLEEP_ASLEEP;
         i++) {
        (void)sleep_step(&fsm, calm);
    }

    SleepWakeWatch watch;
    sleep_wake_watch_begin(&watch, config);

    (void)sleep_step(&fsm, hit);
    sleep_wake_watch_observe(&watch, fsm.state);

    for (unsigned int i = 0;
         i < sleep_hit_to_awake_bound(config) + 10u &&
         fsm.state != SLEEP_AWAKE;
         i++) {
        (void)sleep_step(&fsm, calm);
        sleep_wake_watch_observe(&watch, fsm.state);
    }

    return sleep_wake_watch_report(&watch, fsm.state);
}

/* ------------------------------------------------------------------ */
/* Self-tests                                                          */
/* ------------------------------------------------------------------ */

static SleepFsm make_state(SleepState state, SleepConfig config) {
    SleepFsm fsm;
    sleep_init(&fsm, config);
    fsm.state = state;
    if (state == SLEEP_DROWSY) fsm.drowsy_left = config.drowsy_frames;
    if (state == SLEEP_WAKING) fsm.wake_left = config.wake_frames;
    return fsm;
}

static void test_independent_entry_edges(void) {
    SleepConfig config = sleep_selected_config();
    const SleepInput calm = {0, 0};
    const SleepInput hit = {1, 0};

    SleepFsm awake = make_state(SLEEP_AWAKE, config);
    check("enter_drowsy_from_awake",
          sleep_step(&awake, calm) == 0 &&
          awake.state == SLEEP_DROWSY &&
          awake.last_edge == SLEEP_EDGE_AWAKE_TO_DROWSY &&
          awake.drowsy_left == config.drowsy_frames);

    SleepFsm drowsy = make_state(SLEEP_DROWSY, config);
    drowsy.drowsy_left = 1u; /* boundary: one calm frame completes the edge */
    check("enter_asleep_from_drowsy",
          sleep_step(&drowsy, calm) == 0 &&
          drowsy.state == SLEEP_ASLEEP &&
          drowsy.last_edge == SLEEP_EDGE_DROWSY_TO_ASLEEP);

    SleepFsm asleep = make_state(SLEEP_ASLEEP, config);
    check("enter_waking_from_asleep",
          sleep_step(&asleep, hit) == 0 &&
          asleep.state == SLEEP_WAKING &&
          asleep.last_edge == SLEEP_EDGE_ASLEEP_TO_WAKING &&
          asleep.wake_left == config.wake_frames);

    SleepFsm waking = make_state(SLEEP_WAKING, config);
    waking.wake_left = 1u; /* boundary: one countdown frame completes it */
    check("enter_awake_from_waking",
          sleep_step(&waking, calm) == 0 &&
          waking.state == SLEEP_AWAKE &&
          waking.last_edge == SLEEP_EDGE_WAKING_TO_AWAKE);
}

static void test_transition_trace(void) {
    SleepConfig config = sleep_selected_config();
    SleepFsm fsm;
    sleep_init(&fsm, config);

    const SleepInput calm = {0, 0};
    const SleepInput active = {0, 1};
    const SleepInput hit = {1, 0};

    check("init_awake", fsm.state == SLEEP_AWAKE);

    /* Enter DROWSY and leave it back to AWAKE on activity. */
    (void)sleep_step(&fsm, calm);
    check("leave_drowsy_to_awake",
          fsm.state == SLEEP_DROWSY &&
          sleep_step(&fsm, active) == 0 &&
          fsm.state == SLEEP_AWAKE &&
          fsm.last_edge == SLEEP_EDGE_DROWSY_TO_AWAKE);

    /* Re-enter DROWSY, then leave it through its timeout into ASLEEP. */
    check("reenter_drowsy_from_awake",
          sleep_step(&fsm, calm) == 0 && fsm.state == SLEEP_DROWSY);

    unsigned int holds = 0;
    while (fsm.state == SLEEP_DROWSY && holds <= config.drowsy_frames + 1u) {
        (void)sleep_step(&fsm, calm);
        holds++;
    }
    check("leave_drowsy_to_asleep",
          fsm.state == SLEEP_ASLEEP &&
          fsm.last_edge == SLEEP_EDGE_DROWSY_TO_ASLEEP &&
          holds == config.drowsy_frames);

    /* Leave ASLEEP by collision and enter WAKING in the same frame. */
    check("leave_asleep_to_waking",
          sleep_step(&fsm, hit) == 0 &&
          fsm.state == SLEEP_WAKING &&
          fsm.last_edge == SLEEP_EDGE_ASLEEP_TO_WAKING);

    /* Leave WAKING after exactly W frames and enter AWAKE. */
    unsigned int wake_frames = 0;
    while (fsm.state == SLEEP_WAKING && wake_frames <= config.wake_frames + 1u) {
        (void)sleep_step(&fsm, calm);
        wake_frames++;
    }
    check("leave_waking_to_awake",
          fsm.state == SLEEP_AWAKE &&
          fsm.last_edge == SLEEP_EDGE_WAKING_TO_AWAKE &&
          wake_frames == config.wake_frames);

    /* Once awake, calm starts the next cycle through DROWSY. */
    check("leave_awake_to_drowsy",
          sleep_step(&fsm, calm) == 0 &&
          fsm.state == SLEEP_DROWSY &&
          fsm.last_edge == SLEEP_EDGE_AWAKE_TO_DROWSY);
}

/* Exhaust all raw 2-bit inputs and all base-4 sequences of length 8. */
static void test_completeness_exhaustive(void) {
    SleepConfig config = {2u, 2u};
    const unsigned int len = 8u;
    const unsigned int states = 4u;
    unsigned long long sequences = 1ull << (2u * len); /* 4^len */
    int all_closed = 1;

    for (unsigned int start = 0; start < states; start++) {
        for (unsigned long long code = 0; code < sequences; code++) {
            SleepFsm fsm = make_state((SleepState)start, config);
            unsigned long long x = code;
            for (unsigned int i = 0; i < len; i++) {
                unsigned int raw = (unsigned int)(x & 3ull);
                x >>= 2;
                SleepInput input = {raw & 1u, (raw >> 1) & 1u};
                if (sleep_step(&fsm, input) != 0 ||
                    (unsigned int)fsm.state > SLEEP_WAKING) {
                    all_closed = 0;
                    break;
                }
            }
            if (!all_closed) break;
        }
        if (!all_closed) break;
    }
    check("all_input_sequences_remain_in_four_states", all_closed);

    SleepFsm corrupted = make_state(SLEEP_AWAKE, config);
    corrupted.state = (SleepState)99;
    check("invalid_encoded_state_recovered",
          sleep_step(&corrupted, (SleepInput){0, 0}) == 0 &&
          corrupted.state == SLEEP_AWAKE &&
          corrupted.last_edge == SLEEP_EDGE_RECOVER_INVALID);

    SleepFsm corrupted_and_bad_config = make_state(SLEEP_AWAKE, config);
    corrupted_and_bad_config.state = (SleepState)-7;
    corrupted_and_bad_config.config.drowsy_frames = 0;
    corrupted_and_bad_config.config.wake_frames = 0;
    check("invalid_state_recovered_even_with_bad_config",
          sleep_step(&corrupted_and_bad_config, (SleepInput){0, 0}) == 0 &&
          (unsigned int)corrupted_and_bad_config.state <= SLEEP_WAKING &&
          corrupted_and_bad_config.state == SLEEP_AWAKE);

    check("all_four_state_names_valid",
          strcmp(sleep_state_name(SLEEP_AWAKE), "AWAKE") == 0 &&
          strcmp(sleep_state_name(SLEEP_DROWSY), "DROWSY") == 0 &&
          strcmp(sleep_state_name(SLEEP_ASLEEP), "ASLEEP") == 0 &&
          strcmp(sleep_state_name(SLEEP_WAKING), "WAKING") == 0);
}

static void test_wake_bound(void) {
    SleepConfig selected = sleep_selected_config();
    SleepWakeReport r = sleep_run_canonical_wake(selected);

    if (!r.within_bound) {
        printf("wake_bound_violation actual_frames=%u bound_frames=%u "
               "stuck_edge=%s\n",
               r.actual_frames, r.bound_frames,
               sleep_edge_name(r.stuck_edge));
    }
    check("selected_bound_is_three_frames",
          sleep_edge_frame_budget(SLEEP_EDGE_ASLEEP_TO_WAKING, selected) == 1u &&
          sleep_edge_frame_budget(SLEEP_EDGE_WAKING_TO_AWAKE, selected) == 2u &&
          sleep_hit_to_awake_bound(selected) == SLEEP_WAKE_BUDGET_FRAMES &&
          r.finished && r.within_bound &&
          r.actual_frames == SLEEP_WAKE_BUDGET_FRAMES &&
          r.stuck_edge == SLEEP_EDGE_NONE);

    /* Main loops keep observing after completion; the first AWAKE must freeze
     * the measured latency and must not become a false bound violation. */
    SleepWakeWatch lingering;
    sleep_wake_watch_begin(&lingering, selected);
    sleep_wake_watch_observe(&lingering, SLEEP_WAKING);
    sleep_wake_watch_observe(&lingering, SLEEP_WAKING);
    sleep_wake_watch_observe(&lingering, SLEEP_AWAKE);
    for (unsigned int extra = 0; extra < 5u; extra++)
        sleep_wake_watch_observe(&lingering, SLEEP_AWAKE);
    sleep_wake_watch_observe(&lingering, SLEEP_DROWSY); /* later calm cycle */
    SleepWakeReport lingering_report =
        sleep_wake_watch_report(&lingering, SLEEP_DROWSY);
    check("wake_monitor_latches_completion_and_ignores_later_frames",
          lingering.completed &&
          lingering_report.finished &&
          lingering_report.within_bound &&
          lingering_report.actual_frames == SLEEP_WAKE_BUDGET_FRAMES &&
          lingering_report.stuck_edge == SLEEP_EDGE_NONE);

    int all_w_match = 1;
    for (unsigned int w = 1; w <= SWEEP_MAX_N; w++) {
        SleepConfig config = {8u, w};
        SleepWakeReport rr = sleep_run_canonical_wake(config);
        if (!rr.within_bound) {
            printf("wake_bound_violation W=%u actual_frames=%u "
                   "bound_frames=%u stuck_edge=%s\n",
                   w, rr.actual_frames, rr.bound_frames,
                   sleep_edge_name(rr.stuck_edge));
        }
        if (!rr.finished || !rr.within_bound ||
            rr.actual_frames != 1u + w ||
            rr.bound_frames != 1u + w) {
            all_w_match = 0;
            break;
        }
    }
    check("measured_wake_latency_equals_graph_bound_for_all_W", all_w_match);

    /* A second collision during WAKING must not restart the ramp. */
    SleepConfig repeat_config = {8u, 4u};
    SleepFsm repeat = make_state(SLEEP_WAKING, repeat_config);
    const SleepInput repeated_hit = {1u, 0u};
    unsigned int repeated_frames = 0;
    while (repeat.state == SLEEP_WAKING && repeated_frames < 10u) {
        (void)sleep_step(&repeat, repeated_hit);
        repeated_frames++;
    }
    check("repeated_hit_in_waking_does_not_extend_bound",
          repeat.state == SLEEP_AWAKE && repeated_frames == 4u);

    /* After the first HIT edge, every possible W-frame input tail still
     * reaches AWAKE in exactly W more frames. */
    int arbitrary_tail_ok = 1;
    SleepConfig arbitrary_config = {8u, 4u};
    unsigned long long tails = 1ull << (2u * arbitrary_config.wake_frames);
    for (unsigned long long code = 0; code < tails; code++) {
        SleepFsm fsm = make_state(SLEEP_ASLEEP, arbitrary_config);
        const SleepInput first_hit = {1u, 0u};
        (void)sleep_step(&fsm, first_hit);

        unsigned long long x = code;
        unsigned int tail_frames = 0;
        while (fsm.state != SLEEP_AWAKE &&
               tail_frames < arbitrary_config.wake_frames) {
            unsigned int raw = (unsigned int)(x & 3ull);
            x >>= 2;
            SleepInput input = {raw & 1u, (raw >> 1) & 1u};
            (void)sleep_step(&fsm, input);
            tail_frames++;
        }
        if (fsm.state != SLEEP_AWAKE ||
            tail_frames != arbitrary_config.wake_frames) {
            arbitrary_tail_ok = 0;
            break;
        }
    }
    check("wake_bound_holds_for_every_input_tail", arbitrary_tail_ok);
}

/* Faulty integrator: collision samples disappear before reaching the FSM. */
static SleepWakeReport run_dispatch_fault(SleepConfig config,
                                          unsigned int dropped_hits) {
    SleepFsm fsm;
    sleep_init(&fsm, config);
    const SleepInput calm = {0, 0};
    const SleepInput hit = {1, 0};

    (void)sleep_step(&fsm, calm);
    for (unsigned int i = 0; i < config.drowsy_frames; i++)
        (void)sleep_step(&fsm, calm);

    SleepWakeWatch watch;
    sleep_wake_watch_begin(&watch, config);

    for (unsigned int frame = 0; frame < 100u && fsm.state != SLEEP_AWAKE;
         frame++) {
        if (dropped_hits > 0u) {
            dropped_hits--;
            /* Simulated wiring/callback fault: the hit never calls step(). */
            fsm.frame++;
            fsm.last_input = SLEEP_INPUT_HIT;
            fsm.last_edge = SLEEP_EDGE_ASLEEP_HOLD;
            sleep_wake_watch_observe(&watch, fsm.state);
        } else {
            (void)sleep_step(&fsm, hit);
            sleep_wake_watch_observe(&watch, fsm.state);
        }
    }
    return sleep_wake_watch_report(&watch, fsm.state);
}

/* Faulty counter persistence makes the WAKING edge finish late. */
static SleepWakeReport run_slow_waking_fault(SleepConfig config,
                                             unsigned int extra_frames) {
    SleepFsm fsm;
    sleep_init(&fsm, config);
    const SleepInput calm = {0, 0};
    const SleepInput hit = {1, 0};

    (void)sleep_step(&fsm, calm);
    for (unsigned int i = 0; i < config.drowsy_frames; i++)
        (void)sleep_step(&fsm, calm);

    SleepWakeWatch watch;
    sleep_wake_watch_begin(&watch, config);
    int fault_added = 0;

    for (unsigned int frame = 0; frame < 100u && fsm.state != SLEEP_AWAKE;
         frame++) {
        (void)sleep_step(&fsm, hit);
        if (!fault_added && fsm.state == SLEEP_WAKING) {
            fsm.wake_left += extra_frames;
            fault_added = 1;
        }
        sleep_wake_watch_observe(&watch, fsm.state);
    }
    return sleep_wake_watch_report(&watch, fsm.state);
}

static void test_violation_diagnostics(void) {
    SleepConfig config = sleep_selected_config();

    SleepWakeReport dropped = run_dispatch_fault(config, 3u);
    if (!dropped.within_bound) {
        printf("wake_bound_violation actual_frames=%u bound_frames=%u "
               "stuck_edge=%s\n",
               dropped.actual_frames, dropped.bound_frames,
               sleep_edge_name(dropped.stuck_edge));
    }
    check("dispatch_fault_reports_actual_delay_and_stuck_edge",
          dropped.finished && !dropped.within_bound &&
          dropped.actual_frames == 6u &&
          dropped.bound_frames == 3u &&
          dropped.stuck_edge == SLEEP_EDGE_ASLEEP_TO_WAKING);

    SleepWakeReport slow = run_slow_waking_fault(config, FAULT_DEMONSTRATION_EXTRA);
    if (!slow.within_bound) {
        printf("wake_bound_violation actual_frames=%u bound_frames=%u "
               "stuck_edge=%s\n",
               slow.actual_frames, slow.bound_frames,
               sleep_edge_name(slow.stuck_edge));
    }
    check("slow_waking_fault_reports_actual_delay_and_stuck_edge",
          slow.finished && !slow.within_bound &&
          slow.actual_frames == 3u + FAULT_DEMONSTRATION_EXTRA &&
          slow.bound_frames == 3u &&
          slow.stuck_edge == SLEEP_EDGE_WAKING_TO_AWAKE);
}

/* ------------------------------------------------------------------ */
/* Threshold sweep                                                    */
/* ------------------------------------------------------------------ */

/*
 * Deterministic replay: transient calm gaps of q frames with the listed
 * frequencies.  Gaps shorter than TRUTH_Q are labelled "not safely asleep":
 * entering ASLEEP before their ending collision counts as a false sleep.
 */
static const struct GapFrequency {
    unsigned int calm_frames;
    unsigned int frequency;
} GAP_FREQUENCY[] = {
    {1u, 8u}, {2u, 10u}, {3u, 12u}, {4u, 10u}, {5u, 8u},
    {6u, 6u}, {7u, 4u},  {8u, 3u},  {9u, 2u}
};

#define TRUTH_CALM_FRAMES 10u

static unsigned int false_sleeps_for_d(unsigned int d_frames) {
    SleepConfig config = {d_frames, 2u};
    unsigned int false_sleeps = 0;
    const SleepInput calm = {0, 0};

    for (size_t g = 0; g < sizeof(GAP_FREQUENCY) / sizeof(GAP_FREQUENCY[0]); g++) {
        unsigned int q = GAP_FREQUENCY[g].calm_frames;
        for (unsigned int repeat = 0; repeat < GAP_FREQUENCY[g].frequency; repeat++) {
            SleepFsm fsm;
            sleep_init(&fsm, config);
            int slept_before_hit = 0;

            for (unsigned int i = 0; i < q; i++) {
                (void)sleep_step(&fsm, calm);
                if (fsm.state == SLEEP_ASLEEP) slept_before_hit = 1;
            }
            if (slept_before_hit && q < TRUTH_CALM_FRAMES)
                false_sleeps++;
        }
    }
    return false_sleeps;
}

static int sweep_csv(const char* path) {
    FILE* out = fopen(path, "w");
    if (out == 0) {
        perror(path);
        return 1;
    }

    fprintf(out,
            "candidate,d_false_sleeps,w_measured_wake_frames,"
            "w_bound_frames,coupled_false_sleeps,coupled_wake_frames,"
            "coupled_bound_frames\n");

    for (unsigned int n = 1; n <= SWEEP_MAX_N; n++) {
        unsigned int d_false = false_sleeps_for_d(n);

        SleepConfig w_config = {8u, n};
        SleepWakeReport wake = sleep_run_canonical_wake(w_config);

        unsigned int coupled_false = false_sleeps_for_d(n);
        unsigned int coupled_wake = 1u + n;

        fprintf(out, "%u,%u,%u,%u,%u,%u,%u\n",
                n,
                d_false,
                wake.actual_frames,
                sleep_hit_to_awake_bound(w_config),
                coupled_false,
                coupled_wake,
                coupled_wake);
    }

    if (fclose(out) != 0) {
        perror(path);
        return 1;
    }
    printf("wrote %s\n", path);
    return 0;
}

int main(int argc, char** argv) {
    if (argc == 3 && strcmp(argv[1], "--sweep") == 0)
        return sweep_csv(argv[2]);

    test_independent_entry_edges();
    test_transition_trace();
    test_completeness_exhaustive();
    test_wake_bound();
    test_violation_diagnostics();

    printf("false_sleep_budget=%u\n", SLEEP_FALSE_SLEEP_BUDGET);
    printf("wake_budget_frames=%u\n", SLEEP_WAKE_BUDGET_FRAMES);
    printf("selected_drowsy_frames=%u\n", SLEEP_SELECTED_DROWSY_FRAMES);
    printf("selected_wake_frames=%u\n", SLEEP_SELECTED_WAKE_FRAMES);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
