#!/usr/bin/env python3
"""Plot the angular-momentum precision decision from am_bench CSV output.

Reads  <art>/am_sensitivity.csv, am_tradeoff.csv, am_design_budget.csv
Writes <art>/am_sensitivity.png  and  <art>/am_tradeoff.png
Style follows the dataviz method: one axis per panel, validated categorical
slots 1-3, status colors paired with marker shape + direct labels.
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

# ---- palette (validated light slots 1-3 + fixed status) ----
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
GOOD, CRIT = "#0ca30c", "#d03b3b"
INK, SEC, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, BASE, SURFACE = "#e1e0d9", "#c3c2b7", "#fcfcfb"

for cand in ("Noto Sans CJK SC", "WenQuanYi Zen Hei", "DejaVu Sans"):
    if any(cand in f.name for f in font_manager.fontManager.ttflist):
        plt.rcParams["font.family"] = cand
        break
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 11


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def style_ax(ax):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(BASE)
    ax.tick_params(colors=SEC, labelsize=10)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)


def default_art_dir():
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.normpath(os.path.join(here, "..", "..", "artifacts"))


def plot_sensitivity(art):
    rows = read_csv(os.path.join(art, "am_sensitivity.csv"))
    h = [float(r["h_mm"]) for r in rows]
    canon = [float(r["canonical_rel"]) for r in rows]
    bound = [float(r["analytic_bound_rel"]) for r in rows]
    pmax = [float(r["pair_max_rel"]) for r in rows]
    rms = [float(r["pair_rms_rel"]) for r in rows]
    m = [i for i, x in enumerate(h) if x > 0]

    fig, ax = plt.subplots(figsize=(8.8, 5.4), dpi=170)
    fig.patch.set_facecolor(SURFACE)
    style_ax(ax)

    ax.axhspan(1e-6, 1e-3, color=GOOD, alpha=0.07, zorder=0)
    ax.axhline(1e-3, color=CRIT, linestyle="--", linewidth=2.0, zorder=2)
    ax.text(0.021, 1.45e-3, "千分之一红线  $10^{-3}$", color=CRIT,
            fontsize=10.5, va="bottom", ha="left", fontweight="bold", zorder=6)

    ax.plot([h[i] for i in m], [bound[i] for i in m], "--", color=SEC,
            linewidth=1.8, zorder=3, label="解析最坏上界  $J\\,\\delta r / L_{ref}$")
    ax.plot([h[i] for i in m], [canon[i] for i in m], "o-", color=BLUE,
            linewidth=2.2, markersize=7, zorder=4,
            label="单臂最坏方向（C 实测，与上界重合）")
    ax.plot([h[i] for i in m], [pmax[i] for i in m], "s-", color=ORANGE,
            linewidth=1.6, markersize=5.5, alpha=0.9, zorder=3,
            label="双臂独立误差·3 万例最大")
    ax.plot([h[i] for i in m], [rms[i] for i in m], "^-", color=AQUA,
            linewidth=1.8, markersize=6, zorder=3,
            label="双臂独立误差·RMS")

    drmax = 0.5025
    ax.axvline(drmax, color=BLUE, linestyle=":", linewidth=1.6, alpha=0.8)
    ax.annotate("反解最大力臂误差\n$\\delta r_{max}\\approx0.50$ mm（典型算例）",
                xy=(drmax, 1e-3), xytext=(0.16, 2.4e-5),
                fontsize=10, color=INK,
                arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.4),
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=BASE, lw=0.8))

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.018, 6)
    ax.set_ylim(8e-6, 0.12)
    ax.set_xlabel("力臂误差幅度  $\\delta r$（mm，对数轴）", color=INK, fontsize=11.5)
    ax.set_ylabel("总角动量相对变化  $|\\Delta L| / L_{ref}$（对数轴）",
                   color=INK, fontsize=11.5)
    ax.set_title("力臂每差 1 mm → 相对角动量误差增加  $J\\times10^{-3}/L_{ref}$"
                 " = 1.99‰（典型算例）", color=INK, fontsize=12.5, pad=12)
    leg = ax.legend(loc="upper left", fontsize=9.3, framealpha=0.92,
                    facecolor="white", edgecolor=BASE)
    for t in leg.get_texts():
        t.set_color(INK)
    out = os.path.join(art, "am_sensitivity.png")
    fig.savefig(out, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    return out


def fmt_rel(v):
    if v >= 0.01:
        return f"{v:.2f}"
    s = f"{v:.1e}"            # e.g. 1.5e-03 -> 1.5e-3
    return s.replace("e-0", "e-").replace("e+0", "e")


def plot_tradeoff(art):
    rows = read_csv(os.path.join(art, "am_tradeoff.csv"))
    x = [int(r["x"]) for r in rows]
    mx = [float(r["max_rel"]) for r in rows]
    ns = [float(r["ns_per_collision"]) for r in rows]
    passed = [int(r["passes_1e3"]) == 1 for r in rows]
    labels = [
        "粒子解\n无力臂叉积", "臂误差\n3 mm", "臂误差\n1 mm",
        "臂误差\n0.5 mm", "臂误差\n0.1 mm", "float32\n全叉积",
        "float64\n全叉积", "float64\n+逐次校验",
    ]

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(11.2, 8.4), dpi=170, sharex=True,
        gridspec_kw={"height_ratios": [1.15, 1.0], "hspace": 0.08})
    fig.patch.set_facecolor(SURFACE)

    for ax in (ax1, ax2):
        style_ax(ax)

    # selected operating band: float64 full + guard
    for ax in (ax1, ax2):
        ax.axvspan(5.55, 7.45, color=GOOD, alpha=0.09, zorder=0)

    # ---- top: angular-momentum deviation curve ----
    ax1.semilogy(x, mx, "-", color=BLUE, linewidth=2.2, zorder=3)
    for xi, v, ok in zip(x, mx, passed):
        ax1.scatter([xi], [v], marker="o" if ok else "X", s=95 if ok else 120,
                    facecolor=GOOD if ok else CRIT, edgecolor=INK,
                    linewidth=0.9, zorder=5)
        ax1.annotate(fmt_rel(v), (xi, v), textcoords="offset points",
                     xytext=(0, 11), ha="center", fontsize=8.8, color=INK)
    ax1.axhline(1e-3, color=CRIT, linestyle="--", linewidth=2.0, zorder=2)
    ax1.text(7.42, 1.12e-3, "红线 $10^{-3}$", color=CRIT, fontsize=10.5,
             va="bottom", ha="right", fontweight="bold")
    ax1.set_ylim(1e-17, 12)
    ax1.set_ylabel("最大相对角动量偏差\n（同一 4096 偏心电池）",
                    color=INK, fontsize=11.5)
    ax1.set_title("用「角动量偏差」与「每次碰撞耗时」两条曲线定力臂精度",
                  color=INK, fontsize=13, pad=12, fontweight="bold")
    ax1.text(6.5, 3.2, "采用：double 全叉积\n＋逐次守恒校验",
             ha="center", va="center", fontsize=10.2, color="#006300",
             fontweight="bold",
             bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=GOOD, lw=1.3))
    ax1.scatter([], [], marker="o", facecolor=GOOD, edgecolor=INK, label="整组达标")
    ax1.scatter([], [], marker="X", facecolor=CRIT, edgecolor=INK, label="超千分之一")
    ax1.plot([], [], color=BLUE, lw=2.2, label="最大相对偏差")
    ax1.legend(loc="lower left", fontsize=9.3, framealpha=0.93,
               facecolor="white", edgecolor=BASE, ncol=3)

    # ---- bottom: time-per-collision curve ----
    ax2.plot(x, ns, "-s", color=ORANGE, linewidth=2.2, markersize=8,
             markerfacecolor=ORANGE, markeredgecolor=INK, markeredgewidth=0.7,
             zorder=4)
    for xi, v in zip(x, ns):
        ax2.annotate(f"{v:.0f} ns", (xi, v), textcoords="offset points",
                     xytext=(0, 10), ha="center", fontsize=9, color=INK)
    ax2.set_ylim(0, 92)
    ax2.set_ylabel("每次碰撞耗时（ns，越小越好）", color=INK, fontsize=11.5)
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, fontsize=9.3, color=INK)
    ax2.set_xlim(-0.55, 7.55)
    ax2.set_xlabel("力臂 / 叉积精度档位（越靠右越精确、越有保障）",
                   color=INK, fontsize=11.5)
    ax2.text(0.012, 0.96,
             "本机实测：单次原始叉积 ≈ 2.6 ns；2D 增加一次力臂叉积边际 ≈ 2–3 ns。\n"
             "只要做了完整叉积，毫米级取精几乎不增耗时（33–36 ns 平线）；\n"
             "再加独立守恒校验约 +34 ns，仍不到 0.1 µs/次。",
             transform=ax2.transAxes, ha="left", va="top", fontsize=9.2,
             color=SEC,
             bbox=dict(boxstyle="round,pad=0.45", fc="white", ec=BASE, lw=0.8))

    out = os.path.join(art, "am_tradeoff.png")
    fig.savefig(out, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    return out


def main():
    art = sys.argv[1] if len(sys.argv) > 1 else default_art_dir()
    if not os.path.isdir(art):
        print(f"art dir not found: {art}", file=sys.stderr)
        return 2
    o1 = plot_sensitivity(art)
    o2 = plot_tradeoff(art)
    print("wrote", o1)
    print("wrote", o2)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
