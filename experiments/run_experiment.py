#!/usr/bin/env python3
"""
Degenerate-normal experiment driver.

Drives the real C library (libdegenerate.so) through ctypes: every normal in
every statistic comes from sphere_contact(), never from a Python reimplementation.

Why an explicit event model
---------------------------
Exact centre coincidence has measure zero under free dynamics, but in a real
engine it happens through spawn-on-spawn, positional-correction collapse and
resting stasis. We therefore inject a Poisson stream of coincidence events
throughout the run, in two flavours:

  head_on  centres coincide, relative velocity is nonzero
           -> the VELOCITY branch resolves the limit, no history needed
  stasis   centres coincide and both bodies are dead still; K zero frames
           (the stasis duration, K ~ Geometric(q)) are written into the
           pair's ring before resolution
           -> recovery needs a window deeper than K, otherwise FIXED_AXIS

Recovery vs H is therefore the empirical tail of the stasis distribution, and
the irreducible FIXED_AXIS events are exactly the pairs whose stasis outlasts
the window (plus the frame-0 cold pair). Three ambient overlap tiers set the
initial geometry; exact-coincidence events fire in every tier.

Outputs JSON under artifacts/data/ and a compact stdout report.
Run: python3 -I run_experiment.py [--frames N] [--quick]
"""
import argparse
import ctypes
import json
import math
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
LIB = os.path.join(ROOT, "degenerate", "libdegenerate.so")
OUT = os.path.join(ROOT, "artifacts", "data")

SOURCE_NAMES = ["position", "velocity", "history", "fixed_axis"]
POSITION, VELOCITY, HISTORY, FIXED = 0, 1, 2, 3
FRAME_BYTES = 6 * 4  # PairFrame: two V3 = six floats

N_BALLS = 24
RADIUS = 0.5
BOX = 10.0
DT = 0.02
GRAVITY = -2.0
RESTITUTION = 0.5

# Coincidence event stream (per run, identical across H for a given tier).
EVENT_PROB = 0.8          # P(an event fires on any given frame)
P_STASIS = 0.6            # share of events that arrive dead still
STASIS_Q = 0.10           # Geometric(q): mean stasis 1/q = 10 frames, long tail
STASIS_MAX = 80           # spans the whole H grid so the tail stays visible
WARMUP = 64               # no injections until every ring could be full
COOLDOWN = 8              # a pair is not reused until its zeros have flushed

# Frame-0 audit pairs: every exact-coincidence case must produce a unit normal.
FRAME0_PAIRS = [
    (0, 1, "head_on"),
    (2, 3, "stasis", 2),   # 2 zero frames already behind it
    (4, 5, "cold"),        # born coincident and still, empty ring
]
TIERS = [
    ("shallow", 0.85),   # ambient centre distance -> overlap depth 0.15
    ("medium", 0.50),    # depth 0.50
    ("deep", 0.15),      # depth 0.85
]
H_GRID = [0, 1, 2, 4, 8, 12, 16, 24, 32, 48, 64]


class V3(ctypes.Structure):
    _fields_ = [("x", ctypes.c_float), ("y", ctypes.c_float),
                ("z", ctypes.c_float)]


class Sphere(ctypes.Structure):
    _fields_ = [("c", V3), ("r", ctypes.c_float)]


class Contact(ctypes.Structure):
    _fields_ = [("normal", V3), ("depth", ctypes.c_float),
                ("hit", ctypes.c_int), ("source", ctypes.c_int),
                ("history_age", ctypes.c_int)]


class PairFrame(ctypes.Structure):
    _fields_ = [("rel_pos", V3), ("rel_vel", V3)]


class PairHistory(ctypes.Structure):
    _fields_ = [("ring", ctypes.POINTER(PairFrame)),
                ("capacity", ctypes.c_int), ("head", ctypes.c_int),
                ("count", ctypes.c_int)]


lib = ctypes.CDLL(LIB)
lib.history_init.argtypes = [ctypes.POINTER(PairHistory),
                             ctypes.POINTER(PairFrame), ctypes.c_int]
lib.history_push.argtypes = [ctypes.POINTER(PairHistory), V3, V3]
lib.sphere_contact.argtypes = [ctypes.POINTER(Sphere), ctypes.POINTER(Sphere),
                               V3, ctypes.POINTER(PairHistory)]
lib.sphere_contact.restype = Contact
lib.apply_impulse.argtypes = [ctypes.POINTER(V3), ctypes.POINTER(V3),
                              V3, ctypes.c_float]


def v3(x, y):
    return V3(float(x), float(y), 0.0)


def is_unit(n, tol=1e-5):
    l = math.hypot(n.x, n.y)
    return math.isfinite(n.x) and math.isfinite(n.y) and abs(l - 1.0) <= tol


class Simulation:
    def __init__(self, tier_name, d0, depth, frames, seed):
        self.tier, self.d0, self.H, self.frames = tier_name, d0, depth, frames
        self.rng = random.Random(seed)
        self.pos = [[0.0, 0.0] for _ in range(N_BALLS)]
        self.vel = [[0.0, 0.0] for _ in range(N_BALLS)]
        self._layout()
        self.hist, self.store = {}, {}
        for i in range(N_BALLS):
            for j in range(i + 1, N_BALLS):
                store = (PairFrame * max(depth, 1))()
                h = PairHistory()
                lib.history_init(ctypes.byref(h), store, depth)
                self.hist[(i, j)] = h
                self.store[(i, j)] = store
        self.events = self._schedule_events()
        self.stats = dict(evals=0, contacts=0, sources=[0, 0, 0, 0],
                          failures=[], fixed_events=[],
                          coincidence_audit=[], ambient_degenerate=0,
                          injected=[])     # one record per injected event
        self._seed_frame0()
        self._audit_frame0()

    def _audit_frame0(self):
        """The engineered frame-0 coincidence pairs: every one must return a
        finite unit normal; failures name the ball pair explicitly."""
        s = self.stats
        for spec in FRAME0_PAIRS:
            i, j, kind = spec[0], spec[1], spec[2]
            h = self.hist[(i, j)]
            sa = Sphere(v3(*self.pos[i]), RADIUS)
            sb = Sphere(v3(*self.pos[j]), RADIUS)
            rv = v3(self.vel[j][0] - self.vel[i][0],
                    self.vel[j][1] - self.vel[i][1])
            c = lib.sphere_contact(ctypes.byref(sa), ctypes.byref(sb),
                                   rv, ctypes.byref(h))
            s["coincidence_audit"].append({
                "frame": 0, "ball_a": i, "ball_b": j, "kind": kind,
                "source": SOURCE_NAMES[c.source],
                "normal": [round(c.normal.x, 6), round(c.normal.y, 6)],
                "history_age": c.history_age, "ok": is_unit(c.normal)})
            if not is_unit(c.normal):
                s["failures"].append({
                    "frame": 0, "ball_a": i, "ball_b": j,
                    "source": SOURCE_NAMES[c.source],
                    "normal": [c.normal.x, c.normal.y, c.normal.z]})

    def _layout(self):
        cols, rows = 6, 4
        ox, oy, idx = 1.0, 1.0, 0
        for ry in range(rows):
            for cx in range(cols):
                jx = self.rng.uniform(-0.05, 0.05) * self.d0
                jy = self.rng.uniform(-0.05, 0.05) * self.d0
                self.pos[idx] = [ox + cx * self.d0 + jx,
                                 oy + ry * self.d0 + jy]
                sp = 1.5
                self.vel[idx] = [self.rng.uniform(-sp, sp),
                                 self.rng.uniform(-sp, sp)]
                idx += 1
        for i, j, _kind, *_ in FRAME0_PAIRS:
            x = 1.0 + (i // 2) * 1.8
            # distinct list objects: a chained assignment would alias the
            # two bodies and let the later write clobber the earlier one
            self.pos[i] = [x, BOX - 1.0]
            self.pos[j] = [x, BOX - 1.0]
            self.vel[i] = [0.0, 0.0]
            self.vel[j] = [0.0, 0.0]

    def _schedule_events(self):
        # Drawn here, against the tier RNG only, so the stream is identical
        # for every H in the tier. Warmup leaves rings full of real motion;
        # the cooldown ensures zero frames from one event flush before the
        # same pair is reused, so irreducible FIXED_AXIS events are only the
        # intentional frame-0 cold pair plus genuine K > H tail events.
        out = []
        last_event = {}
        for frame in range(1, self.frames):
            if frame < WARMUP or self.rng.random() >= EVENT_PROB:
                continue
            eligible = [
                (i, j) for i in range(N_BALLS)
                for j in range(i + 1, N_BALLS)
                if frame - last_event.get((i, j), -10 ** 6) > COOLDOWN]
            if not eligible:
                continue
            i, j = self.rng.choice(eligible)
            last_event[(i, j)] = frame
            if self.rng.random() < P_STASIS:
                # geometric draw, 1..STASIS_MAX
                k = min(int(math.log(1.0 - self.rng.random())
                            / math.log(1.0 - STASIS_Q)) + 1, STASIS_MAX)
                out.append((frame, i, j, "stasis", k))
            else:
                out.append((frame, i, j, "head_on", 0))
        return out

    def _seed_frame0(self):
        for spec in FRAME0_PAIRS:
            i, j, kind = spec[0], spec[1], spec[2]
            if kind == "head_on":
                self.vel[j][0] = -2.0
                self.vel[i][0] = 0.5
            elif kind == "stasis":
                # real converging motion in the older slots, then K resting
                # (zero) frames on top: recovery needs a window deeper than K
                h = self.hist[(i, j)]
                for age in range(spec[3] + 1, spec[3] + 7):
                    lib.history_push(ctypes.byref(h),
                                     v3(0.05 * age, 0.0), v3(-1.0, 0.0))
                for _ in range(spec[3]):
                    lib.history_push(ctypes.byref(h),
                                     V3(0.0, 0.0), V3(0.0, 0.0))
            # cold: ring left empty

    def run(self):
        by_frame = {}
        for ev in self.events:
            by_frame.setdefault(ev[0], []).append(ev)
        for frame in range(self.frames):
            self._collide_ambient(frame)
            for ev in by_frame.get(frame, ()):
                self._resolve_event(frame, ev)
            self._integrate()
        return self.stats

    def _resolve_event(self, frame, ev):
        """Controlled-channel resolution of one injected coincidence event.

        Geometry is set exactly on the pair and the normal is read straight
        from the C library, so event metrics are independent of whatever the
        ambient simulation did to these two balls beforehand.
        """
        _f, i, j, kind, k = ev
        h = self.hist[(i, j)]
        if kind == "head_on":
            sa = Sphere(v3(0.0, 0.0), RADIUS)
            sb = Sphere(v3(0.0, 0.0), RADIUS)
            rv = v3(-2.5, 0.0)                 # B approaches A
            c = lib.sphere_contact(ctypes.byref(sa), ctypes.byref(sb),
                                   rv, ctypes.byref(h))
        else:
            for _ in range(k):
                lib.history_push(ctypes.byref(h),
                                 V3(0.0, 0.0), V3(0.0, 0.0))
            sa = Sphere(v3(0.0, 0.0), RADIUS)
            sb = Sphere(v3(0.0, 0.0), RADIUS)
            c = lib.sphere_contact(ctypes.byref(sa), ctypes.byref(sb),
                                   V3(0.0, 0.0), ctypes.byref(h))
        s = self.stats
        s["contacts"] += 1
        s["sources"][c.source] += 1
        if c.source == FIXED:
            s["fixed_events"].append(
                {"frame": frame, "ball_a": i, "ball_b": j, "K": k})
        rec = {"frame": frame, "ball_a": i, "ball_b": j,
               "kind": kind, "K": k, "source": SOURCE_NAMES[c.source],
               "history_age": c.history_age, "ok": is_unit(c.normal)}
        s["injected"].append(rec)
        if not is_unit(c.normal):
            s["failures"].append({
                "frame": frame, "ball_a": i, "ball_b": j,
                "source": SOURCE_NAMES[c.source],
                "normal": [c.normal.x, c.normal.y, c.normal.z]})

    def _collide_ambient(self, frame):
        """Ambient ball-ball contacts: prove the fallback never emits a bad
        normal during ordinary dynamics. Event metrics come from
        _resolve_event instead."""
        s = self.stats
        for i in range(N_BALLS):
            for j in range(i + 1, N_BALLS):
                si = Sphere(v3(*self.pos[i]), RADIUS)
                sj = Sphere(v3(*self.pos[j]), RADIUS)
                rv = v3(self.vel[j][0] - self.vel[i][0],
                        self.vel[j][1] - self.vel[i][1])
                h = self.hist[(i, j)]
                c = lib.sphere_contact(ctypes.byref(si), ctypes.byref(sj),
                                       rv, ctypes.byref(h))
                s["evals"] += 1
                if not c.hit:
                    continue
                if c.source != POSITION:
                    s["ambient_degenerate"] = \
                        s.get("ambient_degenerate", 0) + 1
                if not is_unit(c.normal):
                    s["failures"].append({
                        "frame": frame, "ball_a": i, "ball_b": j,
                        "source": SOURCE_NAMES[c.source],
                        "normal": [c.normal.x, c.normal.y, c.normal.z]})
                    continue
                vi, vj = v3(*self.vel[i]), v3(*self.vel[j])
                lib.apply_impulse(ctypes.byref(vi), ctypes.byref(vj),
                                  c.normal, RESTITUTION)
                self.vel[i] = [vi.x, vi.y]
                self.vel[j] = [vj.x, vj.y]
                if c.depth > 0.0:
                    half = c.depth * 0.5
                    self.pos[i][0] += c.normal.x * half
                    self.pos[i][1] += c.normal.y * half
                    self.pos[j][0] -= c.normal.x * half
                    self.pos[j][1] -= c.normal.y * half

    def _integrate(self):
        for i in range(N_BALLS):
            self.vel[i][1] += GRAVITY * DT
            self.pos[i][0] += self.vel[i][0] * DT
            self.pos[i][1] += self.vel[i][1] * DT
            for axis in (0, 1):
                lo, hi = RADIUS, BOX - RADIUS
                if self.pos[i][axis] < lo:
                    self.pos[i][axis] = lo
                    self.vel[i][axis] = abs(self.vel[i][axis]) * RESTITUTION
                elif self.pos[i][axis] > hi:
                    self.pos[i][axis] = hi
                    self.vel[i][axis] = -abs(self.vel[i][axis]) * RESTITUTION
        for i in range(N_BALLS):
            for j in range(i + 1, N_BALLS):
                lib.history_push(
                    ctypes.byref(self.hist[(i, j)]),
                    v3(self.pos[j][0] - self.pos[i][0],
                       self.pos[j][1] - self.pos[i][1]),
                    v3(self.vel[j][0] - self.vel[i][0],
                       self.vel[j][1] - self.vel[i][1]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frames", type=int, default=500)
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    frames = 120 if args.quick else args.frames
    os.makedirs(OUT, exist_ok=True)

    runs, failures_all, audit_all = [], [], []
    pairs = N_BALLS * (N_BALLS - 1) // 2

    for tier_idx, (tier_name, d0) in enumerate(TIERS):
        for H in H_GRID:
            sim = Simulation(tier_name, d0, H, frames, seed=1000 + tier_idx)
            st = sim.run()
            injected = st["injected"]
            stasis_ev = [e for e in injected if e["kind"] == "stasis"]
            headon_ev = [e for e in injected if e["kind"] == "head_on"]
            fixed_events = [e for e in injected
                            if e["source"] == "fixed_axis"]
            fixed_stasis = [e for e in stasis_ev
                            if e["source"] == "fixed_axis"]
            recovered = [e for e in injected
                         if e["source"] in ("velocity", "history")]
            row = {
                "tier": tier_name, "d0": d0, "H": H,
                "evals": st["evals"], "contacts": st["contacts"],
                "ambient_degenerate": st["ambient_degenerate"],
                "sources": dict(zip(SOURCE_NAMES, st["sources"])),
                "n_injected": len(injected),
                "n_stasis": len(stasis_ev), "n_head_on": len(headon_ev),
                "n_recovered": len(recovered),
                "n_fixed_axis": len(fixed_events),
                # Acceptance metrics over the coincidence event stream.
                # fixed_conditional covers every event; fixed_in_stasis is
                # the one that actually depends on (and bounds) H, because
                # head-on events always resolve from current velocity.
                "fixed_conditional":
                    len(fixed_events) / max(len(injected), 1),
                "fixed_in_stasis":
                    len(fixed_stasis) / max(len(stasis_ev), 1),
                "recovery_rate": len(recovered) / max(len(injected), 1),
                "memory_bytes": pairs * H * FRAME_BYTES,
                "failures": st["failures"],
                "_stasis": stasis_ev if H == max(H_GRID) else [],
            }
            runs.append(row)
            failures_all.extend([dict(f, tier=tier_name, H=H)
                                 for f in st["failures"]])
            if H == max(H_GRID):
                audit_all.extend([dict(a, tier=tier_name)
                                  for a in st["coincidence_audit"]])
            print(f"{tier_name:7s} H={H:3d} ambient_eval={st['evals']:6d} "
                  f"ambient_deg={st['ambient_degenerate']:4d} "
                  f"inj={len(injected):4d} (h={len(headon_ev)},s={len(stasis_ev)}) "
                  f"recovered={len(recovered):4d} fixed={len(fixed_events):4d} "
                  f"fixed_all={row['fixed_conditional']:.4f} "
                  f"fixed_stasis={row['fixed_in_stasis']:.4f} "
                  f"fails={len(st['failures'])}")

    # Stasis duration vs required history age, from the deepest window.
    age_cdf = {}
    for tier_name, _ in TIERS:
        deep = next(r for r in runs if r["tier"] == tier_name
                    and r["H"] == max(H_GRID))
        stasis_ev = deep.pop("_stasis")
        total = len(stasis_ev)
        recovered_ev = [e for e in stasis_ev
                        if e["source"] in ("velocity", "history")]
        # For recovered events history_age should equal K+1; record the
        # fraction of stasis events recoverable with a window of depth H,
        # plus the empirical K distribution we injected.
        k_hist = {k: sum(1 for e in stasis_ev if e["K"] == k)
                  for k in range(1, STASIS_MAX + 1)}
        age_cdf[tier_name] = {
            "samples": total,
            "recovered_at_H64": len(recovered_ev),
            "K_histogram": [{"K": k, "count": k_hist.get(k, 0)}
                            for k in range(1, STASIS_MAX + 1)
                            if k_hist.get(k, 0)],
            "recovery_by_H": [
                {"H": H,
                 "recovered": sum(
                     1 for e in stasis_ev if e["K"] < H),
                 "rate": sum(
                     1 for e in stasis_ev if e["K"] < H) / max(total, 1)}
                for H in H_GRID],
            "age_check": [
                {"K": e["K"], "age": e["history_age"], "source": e["source"]}
                for e in recovered_ev[:25]],
        }
    for r in runs:
        r.pop("_stasis", None)

    # Analytic inversion of the Geometric stasis tail, plus the empirical
    # minimum window satisfying the 1% fixed-axis acceptance line.
    analytic = {
        "q": STASIS_Q, "mean_stasis": 1.0 / STASIS_Q,
        "tail_formula": "P(K >= H) = (1-q)^H, H = history depth",
        "H_for_tail_1pct": math.ceil(math.log(0.01) / math.log(1 - STASIS_Q)),
        "tail_by_H": [{"H": H, "tail": (1 - STASIS_Q) ** H} for H in H_GRID],
    }
    chosen = {}
    for tier_name, _ in TIERS:
        tier_runs = [r for r in runs if r["tier"] == tier_name]
        ok = [r["H"] for r in tier_runs if r["fixed_in_stasis"] <= 0.01]
        chosen[tier_name] = {
            "min_H_for_1pct": min(ok) if ok else None,
            "min_H_for_1pct_all_events": min(
                (r["H"] for r in tier_runs
                 if r["fixed_conditional"] <= 0.01), default=None),
        }
    analytic["empirical_choice"] = chosen

    memory_curve = [
        {"H": H, "bytes_24_balls": pairs * H * FRAME_BYTES,
         "bytes_per_1000_balls": (1000 * 999 // 2) * H * FRAME_BYTES}
        for H in H_GRID]

    with open(os.path.join(OUT, "runs.json"), "w") as f:
        json.dump({"frames": frames, "dt": DT, "n_balls": N_BALLS,
                   "pairs": pairs, "restitution": RESTITUTION,
                   "event_prob": EVENT_PROB, "p_stasis": P_STASIS,
                   "stasis_q": STASIS_Q, "runs": runs}, f, indent=2)
    with open(os.path.join(OUT, "age_cdf.json"), "w") as f:
        json.dump(age_cdf, f, indent=2)
    with open(os.path.join(OUT, "frame0_audit.json"), "w") as f:
        json.dump(audit_all, f, indent=2)
    with open(os.path.join(OUT, "failures.json"), "w") as f:
        json.dump(failures_all, f, indent=2)
    with open(os.path.join(OUT, "memory.json"), "w") as f:
        json.dump({"frame_bytes": FRAME_BYTES, "analytic": analytic,
                   "curve": memory_curve}, f, indent=2)
    print(f"\nwrote JSON to {OUT}")
    print(f"total unusable normals: {len(failures_all)}")


if __name__ == "__main__":
    main()
