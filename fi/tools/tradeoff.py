#!/usr/bin/env python3
"""Coverage-based inversion of the injection-point count.

Two curves decide k, the number of injection points:
  (1) coverage: cumulative functional units covered vs k
  (2) time:     cumulative injection cost vs k
Point k* is where marginal coverage per unit time collapses (the knee), plus the
exact minimum k needed to reach a coverage threshold (solved as an ILP, so the
number is an optimum, not a greedy upper bound).

Universe: every declared functional unit (executable line OR branch arm) in the
modules' functional functions -- collide.test/clampf/closest, solve.resolve,
integrate.step. The modules' own self-tests are excluded; no injection drives them.

Cost model: each injection point is one isolated process (that is how the suite
runs). Measured wall time per process is dominated by a fixed ~0.6 ms exec, so
time is very nearly linear in k. We report both that and the in-process per-call
cost for transparency.

Outputs under results/:
  coverage_matrix.json, coverage_curve.json, time_curve.json, decision.json
  tradeoff.txt (human report)
"""
import json, math, os, subprocess, sys, time, statistics

FI = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOOLS = os.path.join(FI, "tools")
RESULTS = os.path.join(FI, "results")
MODS = ("collide", "solve", "integrate")
sys.path.insert(0, TOOLS)
import paths  # noqa: E402

CORE = [("collide", "N1"), ("collide", "S1"), ("collide", "T1"),
        ("solve", "N2"), ("solve", "S2"), ("solve", "T2"),
        ("integrate", "N3"), ("integrate", "S3"), ("integrate", "T3")]


def covtool(*a):
    return json.loads(subprocess.run([sys.executable, os.path.join(TOOLS, "covtool.py"), *a],
                                     check=True, capture_output=True, text=True).stdout)


def list_cases(m):
    return subprocess.run([os.path.join(FI, "build", m), "list"],
                          check=True, capture_output=True, text=True).stdout.split()


def collect():
    decl, mat = {}, {}
    for m in MODS:
        d = covtool("declared", m)
        decl[m] = {("L", m, x) for x in d["lines"]} | {("B", m, tuple(b)) for b in d["branches"]}
        cases = list_cases(m)
        mx = covtool("matrix", m, *cases)
        for c in cases:
            mat[(m, c)] = ({("L", m, x) for x in mx[c]["lines"]} |
                           {("B", m, tuple(b)) for b in mx[c]["branches"]}) & decl[m]
    return decl, mat


def measure_process_time(m, c, reps=11):
    ts = []
    for _ in range(reps):
        t0 = time.perf_counter()
        subprocess.run([os.path.join(FI, "build", m), c], stdout=subprocess.DEVNULL, check=True)
        ts.append(time.perf_counter() - t0)
    return statistics.median(ts)


def measure_call_time(m, c):
    return covtool("time", m, c, "20000", "60000")["percall_ns"]


def greedy(mat, universe, seed=()):
    order, cov = list(seed), set().union(*[mat[p] for p in seed]) if seed else set()
    rem = set(universe) - cov
    steps = []
    while rem:
        best, gain = None, 0
        for p, u in mat.items():
            if p in order:
                continue
            g = len(u & rem)
            if g > gain:
                gain, best = g, p
        if best is None:
            break
        order.append(best)
        cov |= mat[best]
        rem = set(universe) - cov
        steps.append((best, gain))
    return order, cov, steps


def exact_min_k(mat, universe, threshold):
    """ILP: min #points s.t. >= threshold units covered. Returns (k, chosen)."""
    import numpy as np
    from scipy.optimize import milp, LinearConstraint, Bounds
    U = sorted(universe); idx = {u: i for i, u in enumerate(U)}
    pts = sorted(mat); P, N = len(pts), len(U)
    A = np.zeros((N, P))
    for j, p in enumerate(pts):
        for u in mat[p]:
            A[idx[u], j] = 1
    target = math.ceil(threshold * N)
    c = np.concatenate([np.ones(P), np.zeros(N)])
    M = np.zeros((N, P + N)); M[:, :P] = A
    M[np.arange(N), P + np.arange(N)] = -1
    con1 = LinearConstraint(M, lb=np.zeros(N), ub=np.full(N, np.inf))
    rs = np.zeros((1, P + N)); rs[0, P:] = 1
    con2 = LinearConstraint(rs, lb=np.array([target]), ub=np.array([np.inf]))
    res = milp(c=c, integrality=np.ones(P + N), constraints=[con1, con2],
               bounds=Bounds(np.zeros(P + N), np.ones(P + N)))
    if not res.success:
        raise RuntimeError(res.message)
    x = np.round(res.x[:P]).astype(int)
    chosen = [pts[j] for j in range(P) if x[j]]
    return len(chosen), chosen


def main():
    os.makedirs(RESULTS, exist_ok=True)
    decl, mat = collect()
    universe = set().union(*decl.values())
    N = len(universe)
    assert set().union(*mat.values()) == universe, "universe not fully reachable"

    order, cov, steps = greedy(mat, universe)
    ptime = {p: measure_process_time(*p) for p in mat}
    ctime = {p: measure_call_time(*p) for p in mat}

    # curves over the coverage-first greedy order
    curve, cumt = [], 0.0
    for i, (p, gain) in enumerate(steps, 1):
        cumt += ptime[p]
        c = len(set().union(*[mat[q] for q in order[:i]]))
        curve.append({"k": i, "point": f"{p[0]}/{p[1]}", "gain": gain,
                      "covered": c, "pct": 100.0 * c / N,
                      "cum_ms": cumt * 1e3})

    # knee by max deviation from the chord (0,0)->(k_end, N)
    ke, ye = len(curve), N
    knee = max(curve, key=lambda x: abs(ye * x["k"] - ke * x["covered"]))

    # marginal coverage per added millisecond across the curve
    for i, c in enumerate(curve):
        prev_cov = curve[i - 1]["covered"] if i else 0
        prev_t = curve[i - 1]["cum_ms"] if i else 0.0
        dt = c["cum_ms"] - prev_t
        c["marg_units_per_ms"] = (c["covered"] - prev_cov) / dt if dt > 0 else float("inf")

    # thresholds
    k90, chosen90 = exact_min_k(mat, universe, 0.90)
    k95, chosen95 = exact_min_k(mat, universe, 0.95)
    k100, chosen100 = exact_min_k(mat, universe, 1.00)

    # mandated 3-per-class baseline
    core_cov = set().union(*[mat[p] for p in CORE])
    k_core_aug = None
    _, _, _ = order, cov, steps
    # how many extra greedy points to lift the mandated 9 to 90%
    add, cset = [], set(core_cov)
    rem = universe - cset
    while len(cset) < math.ceil(0.9 * N) and rem:
        best, gain = None, 0
        for p, u in mat.items():
            if p in CORE or p in add:
                continue
            g = len(u & rem)
            if g > gain:
                gain, best = g, p
        if best is None:
            break
        add.append(best); cset |= mat[best]; rem = universe - cset

    decision = {
        "universe_units": N,
        "n_candidate_points": len(mat),
        "knee_k": knee["k"], "knee_pct": knee["pct"],
        "k90": k90, "k95": k95, "k100": k100,
        "mandated_9_units": len(core_cov),
        "mandated_9_pct": 100.0 * len(core_cov) / N,
        "mandated_9_plus": len(add),
        "mandated_9_plus_k": len(CORE) + len(add),
        "mandated_9_plus_pct": 100.0 * len(cset) / N,
        "time_model": "one isolated process per point; ~0.6ms fixed exec",
        "median_process_ms": statistics.median(ptime.values()) * 1e3,
    }

    json.dump({f"{m}/{c}": sorted(map(str, u)) for (m, c), u in mat.items()},
              open(os.path.join(RESULTS, "coverage_matrix.json"), "w"), indent=1)
    json.dump([{"module": p[0], "case": p[1], "process_ms": ptime[p] * 1e3,
                "call_ns": ctime[p]} for p in mat],
              open(os.path.join(RESULTS, "point_times.json"), "w"), indent=1)
    json.dump(curve, open(os.path.join(RESULTS, "coverage_curve.json"), "w"), indent=1)
    json.dump(decision, open(os.path.join(RESULTS, "decision.json"), "w"), indent=1)
    json.dump({"k90": [f"{m}/{c}" for m, c in chosen90],
               "k95": [f"{m}/{c}" for m, c in chosen95],
               "k100": [f"{m}/{c}" for m, c in chosen100]},
              open(os.path.join(RESULTS, "optimal_sets.json"), "w"), indent=1)

    # human-facing path map: which named paths each point reaches
    pathmap = {}
    for m in MODS:
        d = covtool("declared", m)
        cases = list_cases(m)
        mx = covtool("matrix", m, *cases)
        for c in cases:
            pathmap[f"{m}/{c}"] = paths.paths_for_lines(m, mx[c]["lines"])
    json.dump(pathmap, open(os.path.join(RESULTS, "point_paths.json"), "w"), indent=1)

    # path-granularity summary (why it is NOT the 90% metric)
    p_universe = {m: set(paths.all_paths(m, covtool("declared", m)["lines"])) for m in MODS}
    p_all = set().union(*p_universe.values())
    p_core = set().union(*[set(pathmap[f"{m}/{c}"]) & p_universe[m] for m, c in CORE])
    p_curve, p_sel, p_cov = [], [], set()
    while len(p_cov) < len(p_all):
        best, gain = None, 0
        for p, labs in pathmap.items():
            if p in p_sel:
                continue
            g = len(set(labs) - p_cov)
            if g > gain:
                gain, best = g, p
        if best is None:
            break
        p_sel.append(best); p_cov |= set(pathmap[best])
        p_curve.append((len(p_sel), len(p_cov)))
    p90 = next((k for k, c in p_curve if c >= math.ceil(0.9 * len(p_all))), None)
    decision["paths_universe"] = len(p_all)
    decision["paths_core9_pct"] = 100.0 * len(p_core) / len(p_all)
    decision["paths_k_for_90"] = p90
    json.dump(decision, open(os.path.join(RESULTS, "decision.json"), "w"), indent=1)

    L = []
    L.append(f"coverage universe: {N} functional units over {len(mat)} candidate points")
    L.append(f"median per-point process time: {decision['median_process_ms']:.3f} ms")
    L.append("")
    L.append(f"{'k':>3} {'point':>18} {'+units':>7} {'covered':>9} {'pct':>7} "
             f"{'cum_ms':>8} {'u/ms':>8}")
    for c in curve:
        L.append(f"{c['k']:>3} {c['point']:>18} {c['gain']:>7} "
                 f"{c['covered']:>4}/{N:<4} {c['pct']:>6.1f}% {c['cum_ms']:>8.2f} "
                 f"{c['marg_units_per_ms']:>8.1f}")
    L.append("")
    L.append(f"knee of the coverage curve: k={knee['k']} at {knee['pct']:.1f}%")
    L.append(f"exact minimum k >= 90%: {k90}")
    L.append(f"exact minimum k >= 95%: {k95}")
    L.append(f"exact minimum k = 100%: {k100}")
    L.append("")
    L.append(f"mandated 3-per-class (k=9): {len(core_cov)}/{N} = "
             f"{decision['mandated_9_pct']:.1f}%")
    L.append(f"  + {len(add)} augmenting points -> k={len(CORE)+len(add)}, "
             f"{len(cset)}/{N} = {decision['mandated_9_plus_pct']:.1f}%")
    L.append("")
    L.append(f"path granularity (named regions): universe {len(p_all)} paths; "
             f"mandated 9 already reach {decision['paths_core9_pct']:.1f}%; "
             f"k for 90% paths = {p90}. Too coarse to discriminate, so the "
             f"unit metric above is the one k is chosen on.")
    report = "\n".join(L)
    open(os.path.join(RESULTS, "tradeoff.txt"), "w").write(report + "\n")
    print(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
