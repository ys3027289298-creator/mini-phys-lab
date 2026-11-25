#!/usr/bin/env python3
"""Order-sensitive normal/friction impulse experiment.

The model is deliberately the scalar sphere-plane model used by solve/src/solve.c:
unit mass, n perpendicular to t, incoming velocity -u*n + s*t.

Two order twins are compared:

N -> T: snapshot the uncorrected incoming normal load u and use it as the
       friction decrement budget.  This is the conservative formula already
       present in solve.c.

T -> N: before the normal row has produced the actual j_n, use the only
       restitution-independent upper bound, perfect reversal 2*u, as the
       working budget.  If the actual restitution is e < 1, the true normal
       impulse in velocity units is (1+e)*u < 2*u, so the raw tangent pass can
       be outside the Coulomb cone.

A final cone projection gives the admissible result used in the acceptance
table.  All calculations are deterministic.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "friction_order"
OUT.mkdir(parents=True, exist_ok=True)

RUN_MU = 0.30
SWEEP_STEP = 0.01

# Five oblique impacts.  u = -v_n^- > 0, s = |v_t^-| > 0.
# e is deliberately varied so the T-first max-restitution estimate can be
# checked against several true cone widths.
CASES = [
    {"id": 1, "u": 10.0, "s": 8.0, "e": 0.50},
    {"id": 2, "u": 10.0, "s": 7.0, "e": 0.60},
    {"id": 3, "u": 10.0, "s": 6.5, "e": 0.75},
    {"id": 4, "u": 10.0, "s": 6.1, "e": 0.90},
    {"id": 5, "u": 10.0, "s": 3.0, "e": 1.00},
]

BLUE = "#2a78d6"
ORANGE = "#eb6834"
RED = "#d03b3b"
GREEN = "#0ca30c"
GRID = "#e1e0d9"
INK = "#0b0b0b"
MUTED = "#898781"


def clamp(x: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, x))


def accuracy(raw_reduction: float, reference_reduction: float) -> float:
    """100% = exact analytical cone-projected tangential reduction."""
    if reference_reduction <= 1e-12:
        return 100.0
    return 100.0 * clamp(
        1.0 - abs(raw_reduction - reference_reduction) / reference_reduction,
        0.0,
        1.0,
    )


def solve_case(case: dict, mu: float) -> dict:
    u = float(case["u"])
    s = float(case["s"])
    e = float(case["e"])

    normal_correction = (1.0 + e) * u       # j_n / m
    true_cone_budget = mu * normal_correction
    conservative_budget = mu * u
    tangent_first_estimate = 2.0 * mu * u   # max possible reversal, e_max = 1

    reference_dt = min(s, true_cone_budget)

    normal_dt = min(s, conservative_budget)
    tangent_dt_raw = min(s, tangent_first_estimate)

    tangent_dt_projected = min(tangent_dt_raw, true_cone_budget)
    normal_dt_projected = min(normal_dt, true_cone_budget)

    def cone(dt: float) -> dict:
        utilization = dt / true_cone_budget if true_cone_budget > 1e-15 else 0.0
        excess_fraction = max(0.0, utilization - 1.0)
        return {
            "utilization": utilization,
            "excess_pct": 100.0 * excess_fraction,
            "slope": dt / normal_correction if normal_correction > 0 else 0.0,
            "ok": dt <= true_cone_budget + 1e-10,
        }

    n_cone = cone(normal_dt)
    t_cone = cone(tangent_dt_raw)
    np_cone = cone(normal_dt_projected)
    tp_cone = cone(tangent_dt_projected)

    return {
        **case,
        "mu": mu,
        "normal_vn_after": e * u,
        "jn_per_mass": normal_correction,
        "true_cone_budget": true_cone_budget,
        "reference_dt": reference_dt,
        "normal_dt": normal_dt,
        "normal_vt_after": s - normal_dt,
        "normal_acc_pct": accuracy(normal_dt, reference_dt),
        "normal_utilization": n_cone["utilization"],
        "normal_excess_pct": n_cone["excess_pct"],
        "normal_cone_ok": n_cone["ok"],
        "normal_projected_dt": normal_dt_projected,
        "normal_projected_vt_after": s - normal_dt_projected,
        "normal_projected_ok": np_cone["ok"],
        "tangent_dt_raw": tangent_dt_raw,
        "tangent_vt_after_raw": s - tangent_dt_raw,
        "tangent_acc_pct": accuracy(tangent_dt_raw, reference_dt),
        "tangent_utilization": t_cone["utilization"],
        "tangent_excess_pct": t_cone["excess_pct"],
        "tangent_cone_ok_raw": t_cone["ok"],
        "tangent_projected_dt": tangent_dt_projected,
        "tangent_projected_vt_after": s - tangent_dt_projected,
        "tangent_projected_ok": tp_cone["ok"],
    }


def write_five_cases(rows: list[dict]) -> None:
    fields = list(rows[0].keys())
    with (OUT / "five_cases.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    # A compact table is easier to quote in the markdown report.
    with (OUT / "five_cases_table.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "组", "u", "s", "e", "mu",
            "N削减", "N后切向", "N准确率%", "N超出%",
            "T原始削减", "T原始后切向", "T准确率%", "T超出%",
            "T投影后削减", "T投影后切向", "投影后超出%",
        ])
        for r in rows:
            w.writerow([
                r["id"], r["u"], r["s"], r["e"], f'{r["mu"]:.2f}',
                f'{r["normal_dt"]:.4f}', f'{r["normal_vt_after"]:.4f}',
                f'{r["normal_acc_pct"]:.2f}', f'{r["normal_excess_pct"]:.2f}',
                f'{r["tangent_dt_raw"]:.4f}', f'{r["tangent_vt_after_raw"]:.4f}',
                f'{r["tangent_acc_pct"]:.2f}', f'{r["tangent_excess_pct"]:.2f}',
                f'{r["tangent_projected_dt"]:.4f}',
                f'{r["tangent_projected_vt_after"]:.4f}',
                "0.00",
            ])


def write_sweep() -> list[dict]:
    sweep_rows = []
    n_steps = int(round(1.0 / SWEEP_STEP)) + 1
    for k in range(n_steps):
        mu = k * SWEEP_STEP
        run_rows = [solve_case(c, mu) for c in CASES]
        active = mu > 1e-12
        n_violations = sum(active and not r["normal_cone_ok"] for r in run_rows)
        t_violations = sum(active and not r["tangent_cone_ok_raw"] for r in run_rows)
        n_acc = sum(r["normal_acc_pct"] for r in run_rows) / len(run_rows)
        t_acc = sum(r["tangent_acc_pct"] for r in run_rows) / len(run_rows)
        sweep_rows.append({
            "mu": mu,
            "normal_first_violations": n_violations,
            "tangent_first_violations": t_violations,
            "normal_first_accuracy_pct": n_acc,
            "tangent_first_accuracy_pct": t_acc,
            "projected_violations": 0,
        })

    with (OUT / "order_sweep.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(sweep_rows[0].keys()))
        w.writeheader()
        w.writerows(sweep_rows)
    return sweep_rows


def make_plot(sweep: list[dict], five: list[dict]) -> None:
    plt.rcParams.update({
        "font.family": ["Noto Sans CJK SC", "WenQuanYi Zen Hei", "DejaVu Sans"],
        "axes.unicode_minus": False,
        "font.size": 11,
        "axes.edgecolor": "#c3c2b7",
        "axes.labelcolor": INK,
        "text.color": INK,
        "xtick.color": "#52514e",
        "ytick.color": "#52514e",
        "figure.facecolor": "#fcfcfb",
        "axes.facecolor": "#fcfcfb",
    })

    mu = [r["mu"] for r in sweep]
    nv = [r["normal_first_violations"] for r in sweep]
    tv = [r["tangent_first_violations"] for r in sweep]
    na = [r["normal_first_accuracy_pct"] for r in sweep]
    ta = [r["tangent_first_accuracy_pct"] for r in sweep]

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(10.5, 8.2), sharex=True,
        gridspec_kw={"height_ratios": [1, 1.25]},
    )
    fig.subplots_adjust(left=0.09, right=0.97, top=0.91, bottom=0.13, hspace=0.14)
    fig.suptitle("法向优先 vs 切向优先：摩擦锥安全次数与减速准确度", fontsize=15, fontweight="bold", y=0.965)

    for ax in (ax1, ax2):
        ax.grid(True, axis="y", color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

    ax1.plot(mu, nv, color=BLUE, linewidth=2.6, label="法向优先：原始超限数", drawstyle="steps-post")
    ax1.plot(mu, tv, color=ORANGE, linewidth=2.6, label="切向优先：原始超限数", drawstyle="steps-post")
    ax1.axvline(RUN_MU, color=MUTED, linewidth=1.4, linestyle="--")
    run_t_violations = sum(not r["tangent_cone_ok_raw"] for r in five)
    ax1.scatter([RUN_MU], [run_t_violations], s=80, color=RED, zorder=5, edgecolor="white", linewidth=1.5)
    ax1.annotate(
        f"μ={RUN_MU:.2f}\nT-first {run_t_violations}/5 超出",
        xy=(RUN_MU, run_t_violations), xytext=(RUN_MU + 0.08, 3.35),
        arrowprops={"arrowstyle": "->", "color": RED, "lw": 1.4},
        color="#7a1f1f", fontsize=10.5,
    )
    ax1.set_ylabel("五组中的超出次数")
    ax1.set_ylim(-0.25, 5.25)
    ax1.set_yticks(range(0, 6))
    ax1.legend(loc="upper right", frameon=False, ncol=1)

    ax2.plot(mu, na, color=BLUE, linewidth=2.8, label="法向优先平均减速准确度")
    ax2.plot(mu, ta, color=ORANGE, linewidth=2.8, label="切向优先平均减速准确度")
    ax2.axvline(RUN_MU, color=MUTED, linewidth=1.4, linestyle="--")
    n_acc = sum(r["normal_acc_pct"] for r in five) / len(five)
    t_acc = sum(r["tangent_acc_pct"] for r in five) / len(five)
    ax2.scatter([RUN_MU], [n_acc], s=70, color=BLUE, zorder=5, edgecolor="white", linewidth=1.4)
    ax2.scatter([RUN_MU], [t_acc], s=70, color=ORANGE, zorder=5, edgecolor="white", linewidth=1.4)
    ax2.annotate(f"{n_acc:.1f}%", (RUN_MU, n_acc), textcoords="offset points", xytext=(9, 8), color=BLUE, fontweight="bold")
    ax2.annotate(f"{t_acc:.1f}%", (RUN_MU, t_acc), textcoords="offset points", xytext=(9, -14), color=ORANGE, fontweight="bold")
    ax2.set_xlabel("摩擦系数 μ")
    ax2.set_ylabel("相对锥投影参考解的准确度 (%)")
    ax2.set_xlim(0, 1.0)
    ax2.set_ylim(0, 105)
    ax2.legend(loc="lower right", frameon=False)

    fig.text(
        0.09, 0.025,
        "曲线由固定五组斜撞在 μ∈[0,1] 上扫描得到；硬约束是摩擦锥超限次数，其次比较减速准确度。",
        color="#52514e", fontsize=9.5,
    )
    png = OUT / "friction_order_curves.png"
    fig.savefig(png, dpi=180, facecolor=fig.get_facecolor())
    plt.close(fig)


def write_summary(five: list[dict], sweep: list[dict]) -> dict:
    violating = [
        {
            "case": r["id"],
            "excess_pct": round(r["tangent_excess_pct"], 2),
            "utilization": round(r["tangent_utilization"], 6),
        }
        for r in five if not r["tangent_cone_ok_raw"]
    ]
    summary = {
        "run_mu": RUN_MU,
        "case_count": len(five),
        "normal_first_violations": sum(not r["normal_cone_ok"] for r in five),
        "tangent_first_raw_violations": len(violating),
        "tangent_first_violating_cases": violating,
        "all_projected_results_in_cone": all(
            r["normal_projected_ok"] and r["tangent_projected_ok"] for r in five
        ),
        "mean_accuracy_normal_pct": sum(r["normal_acc_pct"] for r in five) / len(five),
        "mean_accuracy_tangent_raw_pct": sum(r["tangent_acc_pct"] for r in five) / len(five),
        "sweep_points": len(sweep),
        "verdict": "normal_first_then_project_tangent_to_mu_jn",
    }
    (OUT / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    five = [solve_case(c, RUN_MU) for c in CASES]
    write_five_cases(five)
    sweep = write_sweep()
    make_plot(sweep, five)
    summary = write_summary(five, sweep)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"artifacts={OUT}")


if __name__ == "__main__":
    main()
