"""
plots.py — 读 results/report.json，产出四类图（PNG）:
  fig1_conflict_curve.png  冲突次数 vs 误判预算：独立 vs 联合（决策曲线 1）
  fig2_tuning_cost.png     调参耗时/评估次数 vs 方案（决策曲线 2）
  fig3_threshold_rates.png 三阈值的解析 vs 经验误判率单调曲线
  fig4_coupling_breakdown.png  C1 推穿 / C2 抖动 / C3 永不休眠 分解
配色取 dataviz 参考调色板前三槽（蓝/橙/青，全配对门限通过）。
"""
import json
import os
import site
import sys

# 绘图脚本不读取任何不可信输入；matplotlib 装在 user site，-I 模式下显式补回
_user_site = site.getusersitepackages()
if _user_site and _user_site not in sys.path:
    sys.path.append(_user_site)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "..", "results")

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
RED, VIOLET = "#e34948", "#6250d6"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e1e0d9"
SURFACE = "#fcfcfb"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Noto Sans CJK SC", "Liberation Sans"],
    "font.size": 10,
    "text.color": INK,
    "axes.edgecolor": "#c3c2b7",
    "axes.labelcolor": INK,
    "axes.titlecolor": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.8,
    "axes.facecolor": SURFACE,
    "figure.facecolor": SURFACE,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def load():
    with open(os.path.join(RESULTS, "report.json")) as f:
        return json.load(f)


def fig1_conflict(report):
    combos = report["combos"]
    alphas = [c["alpha"] for c in combos]
    ind = [c["independent"]["conflicts"] for c in combos]
    joint = [c["joint"]["conflicts"] for c in combos]
    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=150)
    ax.plot(alphas, ind, "-o", color=BLUE, lw=2, ms=8,
            label="independent tuning (closed-form, no coupling term)")
    ax.plot(alphas, joint, "-s", color=ORANGE, lw=2, ms=8,
            label="joint tuning (3-variable constrained search)")
    for x, yi, yj in zip(alphas, ind, joint):
        ax.annotate(f"{yi}", (x, yi), textcoords="offset points",
                    xytext=(0, 9), ha="center", fontsize=9, color=BLUE)
        ax.annotate(f"{yj}", (x, yj), textcoords="offset points",
                    xytext=(0, -14), ha="center", fontsize=9, color=ORANGE)
    ax.set_xlabel("acceptable misjudgement rate  alpha")
    ax.set_ylabel("coupling conflicts per suite  (C1+C2+C3)")
    ax.set_title("Conflict count vs error budget — joint tuning buys fewer "
                 "conflicts at the same alpha", fontsize=11)
    ax.legend(frameon=False, loc="upper right")
    ax.invert_xaxis()
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, "fig1_conflict_curve.png"))
    plt.close(fig)


def fig2_cost(report):
    combos = report["combos"]
    alphas = [c["alpha"] for c in combos]
    t_ind = [c["independent"]["tune_seconds"] for c in combos]
    t_joint = [c["joint"]["tune_seconds"] for c in combos]
    e_joint = [c["joint"]["suite_evals"] for c in combos]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.4), dpi=150)
    width = 0.35
    x = range(len(alphas))
    ax1.bar([i - width / 2 for i in x], t_ind, width, color=BLUE,
            label="independent")
    ax1.bar([i + width / 2 for i in x], t_joint, width, color=ORANGE,
            label="joint")
    ax1.set_xticks(list(x), [f"alpha={a}" for a in alphas])
    ax1.set_ylabel("wall-clock tuning time (s)")
    ax1.set_title("Tuning time: closed-form is instant, joint search ~20 s",
                  fontsize=10.5)
    for i, (a, b) in enumerate(zip(t_ind, t_joint)):
        ax1.annotate(f"{a:.2f}", (i - width / 2, a),
                     textcoords="offset points", xytext=(0, 3),
                     ha="center", fontsize=8.5, color=BLUE)
        ax1.annotate(f"{b:.1f}", (i + width / 2, b),
                     textcoords="offset points", xytext=(0, 3),
                     ha="center", fontsize=8.5, color=ORANGE)
    ax1.legend(frameon=False)

    ax2.plot(alphas, e_joint, "-^", color=AQUA, lw=2, ms=9,
             label="joint: full-suite evaluations")
    ax2.axhline(1, color=BLUE, lw=2, ls="--",
                label="independent: 1 evaluation (analytic)")
    for xv, yv in zip(alphas, e_joint):
        ax2.annotate(str(yv), (xv, yv), textcoords="offset points",
                     xytext=(0, 9), ha="center", fontsize=9, color=AQUA)
    ax2.set_xlabel("acceptable misjudgement rate  alpha")
    ax2.set_ylabel("physics-suite evaluations")
    ax2.set_title("Evaluation cost: 55-60x more simulator runs",
                  fontsize=10.5)
    ax2.invert_xaxis()
    ax2.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, "fig2_tuning_cost.png"))
    plt.close(fig)


def fig3_rates(report):
    titles = {
        "v_sleep": ("sleep threshold v_sleep (m/s)", "false-sleep rate"),
        "d_deg": ("degrade threshold d_deg (m)", "degrade-miss rate"),
        "p_tun": ("penetrate threshold p_tun (m)", "penetration-miss rate"),
    }
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0), dpi=150)
    for ax, which in zip(axes, ("v_sleep", "d_deg", "p_tun")):
        cal = report["calibration"][which]
        ana = report["analytic_curves"][which]
        vals = [r["value"] for r in cal]
        ax.plot(vals, [r["miss_rate"] for r in cal], "-o", color=BLUE,
                lw=2, ms=4, label="empirical (Monte-Carlo calibration)")
        ax.plot(vals, [r["miss"] for r in ana], "--", color=ORANGE,
                lw=2, label="closed-form CDF")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(1e-5, 1.0)
        xlabel, ylabel = titles[which]
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(which, fontsize=10.5)
        ax.legend(frameon=False, fontsize=8.5, loc="lower right")
    fig.suptitle("Monotone threshold-rate relation: empirical inversion "
                 "tracks the analytic CDF (only MC noise separates them)",
                 fontsize=11.5, y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, "fig3_threshold_rates.png"),
                bbox_inches="tight")
    plt.close(fig)


def fig4_coupling(report):
    combos = report["combos"]
    alphas = [c["alpha"] for c in combos]
    parts = ("C1_push_through", "C2_jitter_alarms", "C3_never_slept")
    colors = {"C1_push_through": RED, "C2_jitter_alarms": ORANGE,
              "C3_never_slept": VIOLET}
    labels = {"C1_push_through": "C1 push-through (sleeper driven into slab)",
              "C2_jitter_alarms": "C2 jitter CCD false alarms",
              "C3_never_slept": "C3 slow ball never slept"}

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6), dpi=150,
                             sharey=True)
    for ax, method, title in [
            (axes[0], "independent", "independent tuning"),
            (axes[1], "joint", "joint tuning")]:
        bottoms = [0] * len(combos)
        for p in parts:
            vals = [c[method]["conflict_breakdown"][p] for c in combos]
            ax.bar(range(len(combos)), vals, 0.6, bottom=bottoms,
                   color=colors[p], label=labels[p],
                   edgecolor=SURFACE, linewidth=2)
            for i, (v, btm) in enumerate(zip(vals, bottoms)):
                if v > 0:
                    ax.annotate(str(v), (i, btm + v / 2), ha="center",
                                va="center", fontsize=8.5, color="white",
                                fontweight="bold")
            bottoms = [b + v for b, v in zip(bottoms, vals)]
        ax.set_xticks(range(len(combos)),
                      [f"alpha={a}" for a in alphas])
        ax.set_title(title, fontsize=10.5)
        ax.set_xlabel("")
    axes[0].set_ylabel("conflict events per suite")
    axes[1].legend(frameon=False, fontsize=8.5, loc="upper right")
    fig.suptitle("Coupling conflict decomposition (C1/C2/C3): joint tuning "
                 "removes mostly C2 jitter and C3 never-slept events",
                 fontsize=11.5, y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, "fig4_coupling_breakdown.png"),
                bbox_inches="tight")
    plt.close(fig)


def main():
    report = load()
    fig1_conflict(report)
    fig2_cost(report)
    fig3_rates(report)
    fig4_coupling(report)
    print("figures written to", RESULTS)


if __name__ == "__main__":
    main()
