#!/usr/bin/env python3
"""Render the inclined-ramp figures from the C module's CSV output.

Data is produced by ./out/inclined {curves,tiers,timeseries,sweep}.
Palette: dataviz reference instance, categorical slots 1-3, light surface.
"""
import csv
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogLocator

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "out")
FIGS = os.path.join(HERE, "figs")
os.makedirs(FIGS, exist_ok=True)

# --- palette (validated: all-pairs CVD dE 9.2 light, normal 24.0) -----------
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"   # blue, orange, aqua
REGIME_COLOR = {"rest": S1, "rolling": S2, "sliding": S3}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Noto Sans CJK SC", "DejaVu Sans"],
    "font.size": 10,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "axes.edgecolor": AXIS,
    "axes.labelcolor": INK2,
    "axes.titlecolor": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "text.color": INK,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "lines.linewidth": 2.0,
    "figure.dpi": 200,
})


def read_csv(name):
    with open(os.path.join(OUT, name), newline="") as f:
        return list(csv.DictReader(f))


def style(ax):
    ax.tick_params(length=0)
    ax.set_axisbelow(True)


def save(fig, name):
    path = os.path.join(FIGS, name)
    fig.savefig(path, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    print("wrote", os.path.relpath(path, HERE))


# ---------------------------------------------------------------- fig 1
def fig_critical_angles():
    """How the moment of inertia sets the gap between the two critical angles."""
    mu = 0.389640736
    bodies = [("实心球\nβ=2/5", 0.4), ("空心球\nβ=2/3", 2 / 3), ("圆环\nβ=1", 1.0)]
    ts = math.degrees(math.atan(mu))
    tr = [math.degrees(math.atan(mu * (1 + b) / b)) for _, b in bodies]
    labels = [n for n, _ in bodies]

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    x = range(len(bodies))
    w = 0.34
    b1 = ax.bar([i - w / 2 for i in x], [ts] * len(bodies), w,
                color=S1, label="静止临界角 θ_s", zorder=3)
    b2 = ax.bar([i + w / 2 for i in x], tr, w,
                color=S2, label="滚动临界角 θ_r", zorder=3)
    for bars in (b1, b2):
        for r in bars:
            ax.annotate(f"{r.get_height():.1f}°",
                        (r.get_x() + r.get_width() / 2, r.get_height()),
                        textcoords="offset points", xytext=(0, 3),
                        ha="center", va="bottom", fontsize=9, color=INK2)
    for i, (_, beta) in enumerate(bodies):
        ax.annotate(f"比值 {(1 + beta) / beta:.1f}×",
                    (i, tr[i] + 6.5), ha="center", fontsize=9, color=MUTED)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, color=INK2)
    ax.set_ylim(0, 68)
    ax.set_ylabel("倾角 [°]")
    ax.set_title("两临界角之差由转动惯量决定（μ = 0.390）", pad=12, loc="left")
    ax.legend(frameon=False, loc="upper right", labelcolor=INK2)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.0f}"))
    style(ax)
    save(fig, "fig1_critical_angles.png")


# ---------------------------------------------------------------- fig 2
def fig_calibration():
    rows = read_csv("curves.csv")
    m = [float(r["m"]) for r in rows]
    dev = [float(r["deviation_deg"]) for r in rows]
    fs = [max(float(r["false_slide_rate"]), 1e-6) for r in rows]
    m_star = 1.070529123
    target = 0.005

    fig, (axt, axb) = plt.subplots(2, 1, figsize=(6.6, 5.2), sharex=True,
                                   gridspec_kw={"height_ratios": [1, 1.15]})

    axt.plot(m, dev, color=S1, zorder=3)
    axt.axvline(m_star, color=MUTED, lw=1.2, ls="--", zorder=2)
    axt.plot([m_star], [dev[14]], "o", ms=8, color=S1,
             markeredgecolor=SURFACE, markeredgewidth=2, zorder=4)
    axt.annotate(f"选定 m* = {m_star:.3f}\n偏差 {dev[14]:.2f}°",
                 (m_star, dev[14]), textcoords="offset points", xytext=(12, -26),
                 fontsize=9, color=INK2)
    axt.set_ylabel("临界角偏差 [°]")
    axt.set_title("系数取舍：偏差与滑动误判率", pad=12, loc="left")
    axt.set_ylim(0, max(dev) * 1.10)
    style(axt)

    axb.semilogy(m, fs, color=S2, zorder=3)
    axb.axhline(target, color=S3, lw=1.2, ls=":", zorder=2)
    axb.axvline(m_star, color=MUTED, lw=1.2, ls="--", zorder=2)
    axb.plot([m_star], [target], "o", ms=8, color=S2,
             markeredgecolor=SURFACE, markeredgewidth=2, zorder=4)
    axb.annotate("目标 0.5%", (m[0], target), textcoords="offset points",
                 xytext=(3, 6), fontsize=9, color=S3)
    axb.annotate("m=1.00 时误判 50%\n（恰好压在临界角上）", (m[0], fs[0]),
                 textcoords="offset points", xytext=(8, -4),
                 fontsize=9, color=INK2, va="top")
    axb.annotate("低于 1e-6", (m[-1], fs[-1]), textcoords="offset points",
                 xytext=(-4, 8), ha="right", fontsize=9, color=MUTED)
    axb.set_ylabel("滑动误判率")
    axb.set_xlabel("静摩擦系数倍率 m = μ / tan20°")
    axb.set_xlim(m[0], m[-1])
    axb.set_ylim(1e-6, 1.0)
    axb.yaxis.set_major_locator(LogLocator(base=10, numticks=7))
    axb.yaxis.set_minor_formatter(FuncFormatter(lambda v, _: ""))
    style(axb)

    save(fig, "fig2_calibration.png")


# ---------------------------------------------------------------- fig 3
def fig_regime_sweep():
    rows = read_csv("sweep.csv")
    th = [float(r["theta_deg"]) for r in rows]
    disp = [float(r["disp_m"]) for r in rows]
    slip = [max(float(r["slip_m"]), 1e-10) for r in rows]
    reg = [r["regime"] for r in rows]
    ts, tr = 21.287915, 53.748290
    gate = 0.01

    fig, (axt, axb) = plt.subplots(2, 1, figsize=(7.0, 5.6), sharex=True,
                                   gridspec_kw={"height_ratios": [1, 1.1]})

    for ax in (axt, axb):
        ax.axvspan(0, ts, color=S1, alpha=0.07, lw=0, zorder=0)
        ax.axvspan(ts, tr, color=S2, alpha=0.07, lw=0, zorder=0)
        ax.axvspan(tr, 70, color=S3, alpha=0.07, lw=0, zorder=0)
        ax.axvline(ts, color=INK2, lw=1.0, ls="--", zorder=2)
        ax.axvline(tr, color=INK2, lw=1.0, ls="--", zorder=2)
        style(ax)

    axt.plot(th, disp, color=INK2, zorder=3)
    axt.set_ylabel("10 s 内位移 [m]")
    axt.set_title("倾角扫描：静止 / 滚动 / 滑动三档", pad=12, loc="left")
    axt.text(ts * 0.45, 8, "静止", ha="center", fontsize=10, color=S1)
    axt.text((ts + tr) / 2, 8, "滚动（接触点不滑）", ha="center",
             fontsize=10, color=S2)
    axt.text((tr + 70) / 2, 8, "滑动", ha="center", fontsize=10, color=S3)
    axt.text(ts - 0.8, axt.get_ylim()[1] * 0.90, f"θ_s = {ts:.1f}°",
             ha="right", fontsize=9, color=INK2)
    axt.text(tr - 0.8, axt.get_ylim()[1] * 0.90, f"θ_r = {tr:.1f}°",
             ha="right", fontsize=9, color=INK2)
    i_jump = min(range(len(th)), key=lambda i: abs(th[i] - ts - 1.0))
    axt.annotate("越过 θ_s 即起滚：位移从 0 跳到百米级",
                 (th[i_jump], disp[i_jump]), textcoords="offset points",
                 xytext=(14, 18), fontsize=9, color=INK2)

    axb.semilogy(th, slip, color=INK2, zorder=3)
    axb.axhline(gate, color=S3, lw=1.2, ls=":", zorder=2)
    axb.annotate("验收门 R/100", (0.5, gate), textcoords="offset points",
                 xytext=(0, 5), fontsize=9, color=S3)
    axb.text(5, 1e-9, "θ < θ_r：接触点滑移恰为 0", fontsize=9, color=INK2)
    axb.set_ylabel("10 s 内接触点滑移 [m]")
    axb.set_xlabel("倾角 θ [°]")
    axb.set_xlim(0, 70)
    axb.set_ylim(1e-10, 1e3)
    save(fig, "fig3_regime_sweep.png")


# ---------------------------------------------------------------- fig 4
def fig_tiers():
    rows = read_csv("timeseries.csv")
    by = {}
    for r in rows:
        by.setdefault(r["tier"], []).append((float(r["t"]), float(r["s_m"])))
    order = [("rest", "静止 20.3°", S1),
             ("rolling", "滚动 22.3°", S2),
             ("sliding", "滑动 54.7°", S3)]

    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    for key, label, color in order:
        ts = [p[0] for p in by[key]]
        ss = [p[1] for p in by[key]]
        ax.plot(ts, ss, color=color, label=label, zorder=3)
        yoff = {"rest": 0, "rolling": 0, "sliding": 0}[key]
        ax.annotate(label, (ts[-1], ss[-1]), textcoords="offset points",
                    xytext=(8, yoff), va="center", fontsize=9, color=color)
    ax.set_xlabel("时间 [s]")
    ax.set_ylabel("沿坡位移 s [m]")
    ax.set_title("三档倾角各跑 10 s（R = 1 m，验收门 R/100 = 0.01 m）",
                 pad=12, loc="left")
    ax.set_xlim(0, 11.6)
    ax.set_ylim(-12, 320)
    ax.legend(frameon=False, loc="upper left", labelcolor=INK2)
    style(ax)
    save(fig, "fig4_tiers.png")


if __name__ == "__main__":
    fig_critical_angles()
    fig_calibration()
    fig_regime_sweep()
    fig_tiers()
