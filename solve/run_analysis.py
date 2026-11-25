#!/usr/bin/env python3
"""
run_analysis.py -- reproduce every measured number behind the report.

Runs ./bench_solve modes with a pinned CPU when possible, collects real timer
output into <out>/ (default ../artifacts):

  solve_raw_output.txt   human output of scenario/baseline/final
  baseline.json          three-way phase timing
  final.json             final pass/fail verdict numbers
  iter_time.csv          fixed K vs solve time
  convergence.csv        per-iteration impulse delta (cold stress frame)
  eps_sweep.csv          threshold -> residual / time / iterations
  results.json           everything + least-squares fits + environment

Plotting is separate: make_curves.py
"""
import argparse, csv, json, os, platform, re, subprocess, sys, datetime

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.join(HERE, "bench_solve")

EPS_GRID = [1.0, 0.5, 0.3, 0.2, 0.1, 0.06, 0.04, 0.03, 0.025, 0.02,
            0.015, 0.012, 0.01, 0.008, 0.005, 0.003, 0.001]
FRAMES = 2400


def run(args):
    p = subprocess.run([BENCH] + [str(a) for a in args],
                       capture_output=True, text=True)
    if p.returncode not in (0, 1):  # final returns 1 on failed verdict
        sys.stderr.write(p.stderr)
        raise SystemExit(f"bench_solve {args} failed rc={p.returncode}")
    return p.stdout


def last_json(text):
    m = re.findall(r"^JSON (\{.*\})\s*$", text, re.M)
    if not m:
        raise ValueError("no JSON line in output:\n" + text)
    return json.loads(m[-1])


def fit_linear(xs, ys):
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    b = sxy / sxx
    a = my - b * mx
    ss_res = sum((y - (a + b * x)) ** 2 for x, y in zip(xs, ys))
    ss_tot = sum((y - my) ** 2 for y in ys)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return a, b, r2


def cpu_info():
    info = {"processor": platform.processor(), "python": platform.python_version()}
    try:
        for line in open("/proc/cpuinfo"):
            if line.startswith("model name"):
                info["cpu"] = line.split(":", 1)[1].strip()
                break
    except OSError:
        pass
    return info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "..", "artifacts"))
    ap.add_argument("--baseline-reps", type=int, default=400)
    ap.add_argument("--sweep-reps", type=int, default=150)
    ap.add_argument("--no-pin", action="store_true")
    args = ap.parse_args()
    out = os.path.abspath(args.out)
    os.makedirs(out, exist_ok=True)

    if not os.path.exists(BENCH):
        raise SystemExit("bench_solve missing; run `make bench` first")

    raw = []

    scenario_txt = run(["scenario"])
    raw.append("### bench_solve scenario\n" + scenario_txt)

    base_txt = run(["baseline", 4096, args.baseline_reps])
    raw.append("### bench_solve baseline 4096 %d\n" % args.baseline_reps + base_txt)
    baseline = last_json(base_txt)

    final_txt = run(["final", "0.01", 2001, 1200])
    raw.append("### bench_solve final 0.01 2001 1200\n" + final_txt)
    final = last_json(final_txt)

    iters_txt = run(["itersweep", args.sweep_reps])
    raw.append("### bench_solve itersweep %d\n" % args.sweep_reps + iters_txt)
    iter_path = os.path.join(out, "iter_time.csv")
    with open(iter_path, "w", newline="") as f:
        lines = [ln for ln in iters_txt.splitlines() if ln.strip()]
        f.write("\n".join(lines) + "\n")
    iter_rows = list(csv.DictReader(lines))
    K = [int(r["K"]) for r in iter_rows]
    leg = [float(r["legacy_median_us"]) for r in iter_rows]
    cac = [float(r["cachedfix_median_us"]) for r in iter_rows]
    a_l, b_l, r2_l = fit_linear(K, leg)
    a_c, b_c, r2_c = fit_linear(K, cac)

    conv_txt = run(["converge", 512])
    raw.append("### bench_solve converge 512\n" + conv_txt)
    conv_lines = [ln for ln in conv_txt.splitlines() if ln.strip()]
    with open(os.path.join(out, "convergence.csv"), "w", newline="") as f:
        f.write("\n".join(conv_lines) + "\n")

    # epsilon sweep
    eps_rows = []
    eps_raw = ["eps,res_max,res_steady,res_mean,solve_med_us,solve_p95_us,"
               "it_med,it_max,hit_cap,pass_residual"]
    for eps in EPS_GRID:
        txt = run(["sim", repr(eps), FRAMES])
        d = last_json(txt)
        d["pass_residual"] = int(d["res_max"] < 1e-3)
        eps_rows.append(d)
        eps_raw.append(",".join(str(d[k]) for k in [
            "eps", "res_max", "res_steady", "res_mean", "solve_med_us",
            "solve_p95_us", "it_med", "it_max", "hit_cap", "pass_residual"]))
    with open(os.path.join(out, "eps_sweep.csv"), "w") as f:
        f.write("\n".join(eps_raw) + "\n")
    raw.append("### eps sweep (%d frames each)\n" % FRAMES +
               "\n".join(json.dumps(r) for r in eps_rows))

    # loosest eps that still passes with margin levels
    passing = [r for r in eps_rows if r["pass_residual"]]
    loosest_pass = max(passing, key=lambda r: r["eps"]) if passing else None
    chosen = next((r for r in eps_rows if abs(r["eps"] - 0.01) < 1e-12), None)

    results = {
        "generated_utc": datetime.datetime.utcnow().isoformat() + "Z",
        "environment": cpu_info(),
        "frames_per_eps": FRAMES,
        "baseline": baseline,
        "final": final,
        "fit_legacy_us": {"intercept": a_l, "slope_per_iter": b_l, "r2": r2_l,
                          "K_at_4000us": max(0.0, (4000.0 - a_l) / b_l)},
        "fit_cached_us": {"intercept": a_c, "slope_per_iter": b_c, "r2": r2_c},
        "speedup_outer_vs_final_p50":
            baseline["outer_median_us"] / final["median_us"],
        "speedup_outer_vs_final_p95":
            baseline["outer_p95_us"] / final["p95_us"],
        "loosest_pass_eps": loosest_pass["eps"] if loosest_pass else None,
        "chosen_eps_row": chosen,
        "eps_sweep": eps_rows,
    }
    with open(os.path.join(out, "results.json"), "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    with open(os.path.join(out, "solve_raw_output.txt"), "w") as f:
        f.write("\n\n".join(raw) + "\n")

    print(json.dumps({
        "baseline_outer_median_us": baseline["outer_median_us"],
        "shares": [baseline["enum_share"], baseline["redund_share"],
                   baseline["core_share"]],
        "fit_legacy": [a_l, b_l, r2_l],
        "K_at_4ms": results["fit_legacy_us"]["K_at_4000us"],
        "final_median_us": final["median_us"],
        "final_p95_us": final["p95_us"],
        "res_max": final["res_max"],
        "overall": final["overall"],
        "loosest_pass_eps": results["loosest_pass_eps"],
    }, indent=2))


if __name__ == "__main__":
    main()
