#!/usr/bin/env python3
"""Plot sleep-threshold sweep curves and the three-grade acceptance report.

Usage: plot_sweep.py <sweep.csv> <grade_runs.json> <out_dir>
"""
import csv
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np

# ---- validated reference palette (light surface) -------------------------
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
GRADE_COLORS = {"5": BLUE, "10": ORANGE, "20": AQUA}
INK = "#0b0b0b"
SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"
SURFACE = "#fcfcfb"
SUCCESS_TEXT = "#006300"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Noto Sans CJK SC", "system-ui", "Segoe UI",
                        "DejaVu Sans"],
    "axes.unicode_minus": False,
    "font.size": 10.5,
    "axes.facecolor": SURFACE,
    "figure.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "axes.edgecolor": BASELINE,
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": SECONDARY,
    "ytick.color": SECONDARY,
    "axes.titlecolor": INK,
})


def style_ax(ax):
    ax.grid(True, color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(BASELINE)
        ax.spines[side].set_linewidth(1.0)
    ax.tick_params(length=0)


def load_sweep(path):
    cols = {}
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    for key in rows[0]:
        cols[key] = np.array([float(r[key]) for r in rows])
    return cols


def plot_curves(sw, grade_doc, out_path):
    t_sel = grade_doc["threshold"]
    jitter = grade_doc["jitter_floor"]
    vmin = {str(g["grade_deg"]): g["v_min_window"]
            for g in grade_doc["grades"]}
    t = sw["threshold"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.6, 4.9))

    # safe operating band: above solver jitter, below shallow-grade knee
    for ax in (ax1, ax2):
        ax.axvspan(jitter, vmin["5"], color=AQUA, alpha=0.08, zorder=0)
        ax.axvline(t_sel, color=INK, linestyle=(0, (5, 3)), linewidth=1.8,
                   zorder=3)

    # ---- left: false sleep count -----------------------------------------
    for deg in ("5", "10", "20"):
        y = sw[f"false_{deg}"]
        ax1.plot(t, y, drawstyle="steps-post", color=GRADE_COLORS[deg],
                 linewidth=2.2, marker="o", markersize=3.6,
                 markerfacecolor=SURFACE, markeredgewidth=1.4,
                 markevery=2, label=f"{deg}° 坡", zorder=4)
        # analytic confirmation-window bound
        ax1.axvline(vmin[deg], color=GRADE_COLORS[deg], linestyle=":",
                    linewidth=1.3, alpha=0.55, zorder=2)
        ax1.annotate(f"v_min={vmin[deg]:.2f}",
                     xy=(vmin[deg], 8.15), xycoords=("data", "data"),
                     ha="center", va="bottom", fontsize=8.6,
                     color=GRADE_COLORS[deg])

    ax1.text(t_sel, 3.9, f"选定 T={t_sel:g}", rotation=90,
             ha="center", va="center", fontsize=9.6, fontweight="bold",
             color=INK, zorder=5,
             bbox=dict(boxstyle="round,pad=0.25", fc=SURFACE, ec="none",
                       alpha=0.85))
    ax1.set_xscale("log")
    ax1.set_xlim(t.min() * 0.92, t.max() * 1.08)
    ax1.set_ylim(-0.6, 9.4)
    ax1.set_yticks(range(0, 9, 2))
    ax1.set_xlabel("休眠速度阈值 $v_s$（m/s，对数轴）")
    ax1.set_ylabel("坡上误睡球数（每档 8 个探针）")
    ax1.set_title("① 误判数：阈值越过各坡最小滚动速度后陡增",
                  fontsize=12, pad=10)
    leg1 = ax1.legend(loc="upper left", frameon=False, ncol=3,
                      bbox_to_anchor=(0.0, 1.02))
    for txt in leg1.get_texts():
        txt.set_color(SECONDARY)
    style_ax(ax1)

    # ---- right: sleep ratio + collision skip -----------------------------
    ax2.plot(t, sw["bed_ratio"], color=BLUE, linewidth=2.4, marker="o",
             markersize=4.0, markerfacecolor=SURFACE, markeredgewidth=1.5,
             markevery=2, label="末帧休眠球占比", zorder=4)
    ax2.plot(t, sw["bed_skip"], color=AQUA, linewidth=2.0, linestyle="--",
             label="碰撞对跳过比例（全程）", zorder=3)
    ax2.axvline(jitter, color=MUTED, linestyle=":", linewidth=1.6, zorder=2)
    ax2.annotate(f"抖动下界 {jitter:.3f}", xy=(jitter * 1.04, 0.06),
                 ha="left", va="bottom", fontsize=8.8, color=SECONDARY)
    ax2.annotate(f"选定 T={t_sel:g}", xy=(t_sel, 0.55),
                 xytext=(t_sel * 1.55, 0.42), fontsize=9.5,
                 fontweight="bold", color=INK,
                 arrowprops=dict(arrowstyle="-", color=MUTED, lw=1.0))
    ax2.text(np.sqrt(jitter * vmin["5"]), 0.16, "安全工作区",
             ha="center", va="center", fontsize=9.5, color=SUCCESS_TEXT,
             alpha=0.85)

    ax2.set_xscale("log")
    ax2.set_xlim(t.min() * 0.92, t.max() * 1.08)
    ax2.set_ylim(-0.03, 1.08)
    ax2.set_yticks(np.linspace(0, 1, 6))
    ax2.set_xlabel("休眠速度阈值 $v_s$（m/s，对数轴）")
    ax2.set_ylabel("比例")
    ax2.set_title("② 休眠比例：越过低阈值抖动地板后饱和到 1",
                  fontsize=12, pad=10)
    leg2 = ax2.legend(loc="center right", frameon=False)
    for txt in leg2.get_texts():
        txt.set_color(SECONDARY)
    style_ax(ax2)

    fig.tight_layout(rect=(0, 0, 1, 1))
    fig.savefig(out_path, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path}")


def plot_grades(grade_doc, out_path):
    grades = grade_doc["grades"]
    labels = [f"{int(g['grade_deg'])}° " +
              {5: "缓坡", 10: "中坡", 20: "陡坡"}[int(g["grade_deg"])]
              for g in grades]
    xs = np.arange(len(grades))
    end_r = [g["sleep_ratio"] for g in grades]
    duty = [g["duty_ratio"] for g in grades]
    skip = [g["pair_skip_ratio"] for g in grades]
    colors = [GRADE_COLORS[str(int(g["grade_deg"]))] for g in grades]

    fig, ax = plt.subplots(figsize=(7.6, 4.9))
    w = 0.36
    b1 = ax.bar(xs - w / 2, end_r, w, color=colors, zorder=3,
                label="末帧休眠球占比")
    b2 = ax.bar(xs + w / 2, duty, w, color=colors, alpha=0.30, zorder=3,
                edgecolor=colors, linewidth=1.2, label="平均休眠球帧占比")

    for x, v in zip(xs - w / 2, end_r):
        ax.text(x, v + 0.02, f"{v*100:.1f}%", ha="center", va="bottom",
                fontsize=10.5, fontweight="bold", color=INK)
    for x, v in zip(xs + w / 2, duty):
        ax.text(x, v + 0.02, f"{v*100:.0f}%", ha="center", va="bottom",
                fontsize=9.2, color=SECONDARY)

    for i, g in enumerate(grades):
        assert g["false_sleeps"] == 0 and g["false_oracle"] == 0
        ax.text(xs[i], 1.10, "误睡 0 · oracle 0", ha="center", va="bottom",
                fontsize=9.8, color=SUCCESS_TEXT, fontweight="bold")
        ax.text(xs[i], -0.155, f"碰撞对跳过 {skip[i]*100:.1f}%",
                ha="center", va="top", fontsize=9.0, color=SECONDARY,
                transform=ax.get_xaxis_transform())

    ax.set_xticks(xs)
    ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylim(0, 1.20)
    ax.set_yticks(np.linspace(0, 1, 6))
    ax.set_ylabel("占比")
    ax.set_title(f"选定 T={grade_doc['threshold']:g} m/s、"
                 f"N={grade_doc['confirm_frames']} 帧下三档坡度各跑一次（误判均为 0）",
                 fontsize=12, pad=46)
    leg = ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.005),
                    ncol=2, frameon=False, columnspacing=2.2)
    for txt in leg.get_texts():
        txt.set_color(SECONDARY)
    style_ax(ax)

    fig.tight_layout()
    fig.savefig(out_path, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path}")


def main():
    sweep_path, grades_path, out_dir = sys.argv[1], sys.argv[2], sys.argv[3]
    sw = load_sweep(sweep_path)
    with open(grades_path) as f:
        grade_doc = json.load(f)
    os.makedirs(out_dir, exist_ok=True)
    plot_curves(sw, grade_doc,
                os.path.join(out_dir, "sleep_threshold_curves.png"))
    plot_grades(grade_doc,
                os.path.join(out_dir, "sleep_grade_report.png"))


if __name__ == "__main__":
    main()
