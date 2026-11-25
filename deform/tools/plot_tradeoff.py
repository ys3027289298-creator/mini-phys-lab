#!/usr/bin/env python3
"""
Plot the two stiffness trade-off curves from artifacts/stiffness_sweep.csv:
squash visibility (max squash ratio) and bounce height vs stiffness k.

Two panels sharing one log-x axis — never a dual-axis chart. The calibrated
stiffness is marked on both panels.  Run:

  python3 tools/plot_tradeoff.py artifacts/stiffness_sweep.csv 32792.96 \
      artifacts/stiffness_tradeoff.png
"""
import sys
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# validated dataviz categorical slots 1/2, light-mode chrome
BLUE, ORANGE = "#2a78d6", "#eb6834"
SURFACE, GRID, INK, MUTED = "#fcfcfb", "#e1e0d9", "#0b0b0b", "#898781"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["system-ui", "Segoe UI", "DejaVu Sans"],
    "font.size": 10,
    "text.color": INK,
    "axes.edgecolor": "#c3c2b7",
    "axes.labelcolor": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
})


def load(path):
    ks, eps, bounce = [], [], []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            ks.append(float(row["k_Npm"]))
            eps.append(float(row["squash_ratio"]))
            bounce.append(float(row["bounce_height_m"]))
    return ks, eps, bounce


def pct(x, _):
    return f"{x * 100:.0f}%"


def main():
    csv_path, k_cal, out_path = sys.argv[1], float(sys.argv[2]), sys.argv[3]
    ks, eps, bounce = load(csv_path)

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(8, 7), sharex=True, facecolor=SURFACE)
    for idx, ax in enumerate((ax1, ax2)):
        ax.set_facecolor(SURFACE)
        ax.grid(True, color=GRID, linewidth=0.8)
        ax.axvline(k_cal, color=MUTED, linestyle="--", linewidth=1.2)
        ax.annotate(f"calibrated k = {k_cal:,.0f} N/m",
                    xy=(k_cal, 0.92 if idx == 0 else 0.08),
                    xycoords=("data", "axes fraction"),
                    xytext=(8, 0), textcoords="offset points",
                    color=MUTED, fontsize=9, ha="left",
                    va="top" if idx == 0 else "bottom")

    ax1.plot(ks, eps, color=BLUE, linewidth=2, marker="o",
             markersize=4, label="max squash ratio")
    ax1.set_ylabel("max squash ratio  ε")
    ax1.yaxis.set_major_formatter(FuncFormatter(pct))
    ax1.set_title("Squash visibility vs stiffness — softer = jelly, "
                  "stiffer = invisible", fontsize=11, loc="left")

    ax2.plot(ks, bounce, color=ORANGE, linewidth=2, marker="o",
             markersize=4, label="bounce height")
    ax2.set_ylabel("bounce height [m]")
    ax2.set_xlabel("stiffness k [N/m]")
    ax2.set_title("Bounce height vs stiffness — near-flat at fixed "
                  "restitution; only jelly-soft contacts sag",
                  fontsize=11, loc="left")
    ax2.set_xscale("log")

    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_path, dpi=150, facecolor=SURFACE)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
