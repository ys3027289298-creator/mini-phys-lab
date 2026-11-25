#!/usr/bin/env python3
"""Plot the sleep FSM calibration curves from the C sweep CSV."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
BLUE = "#2a78d6"
ORANGE = "#eb6834"
GREEN = "#1baf7a"

SELECTED_D = 6
SELECTED_W = 2
FALSE_SLEEP_BUDGET = 10
WAKE_BUDGET = 3


def read_rows(path: Path) -> list[dict[str, int]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [{k: int(v) for k, v in row.items()} for row in csv.DictReader(handle)]


def style_axis(ax: plt.Axes) -> None:
    ax.set_facecolor(SURFACE)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=SECONDARY, labelsize=9)
    ax.xaxis.label.set_color(SECONDARY)
    ax.yaxis.label.set_color(SECONDARY)
    ax.title.set_color(INK)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--png", required=True, type=Path)
    args = parser.parse_args()

    rows = read_rows(args.csv)
    n = [r["candidate"] for r in rows]
    false_sleeps = [r["d_false_sleeps"] for r in rows]
    wake = [r["w_measured_wake_frames"] for r in rows]

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [
                "Noto Sans CJK SC",
                "WenQuanYi Zen Hei",
                "system-ui",
                "-apple-system",
                "Segoe UI",
                "sans-serif",
            ],
            "figure.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
        }
    )
    fig, (ax_false, ax_wake) = plt.subplots(1, 2, figsize=(11.2, 4.4), dpi=170)

    # Left: false sleeps versus the drowsy dwell D.
    style_axis(ax_false)
    ax_false.axvspan(SELECTED_D, max(n), color=GREEN, alpha=0.08, linewidth=0)
    ax_false.axhline(
        FALSE_SLEEP_BUDGET,
        color=SECONDARY,
        linestyle=(0, (5, 3)),
        linewidth=1.2,
    )
    ax_false.plot(
        n,
        false_sleeps,
        color=BLUE,
        linewidth=2.2,
        marker="o",
        markersize=5.5,
        markerfacecolor=BLUE,
        markeredgecolor=SURFACE,
        markeredgewidth=1.4,
        label="误睡次数",
    )
    d_value = false_sleeps[n.index(SELECTED_D)]
    ax_false.scatter(
        [SELECTED_D],
        [d_value],
        s=110,
        color=BLUE,
        edgecolor=INK,
        linewidth=1.1,
        zorder=5,
    )
    ax_false.annotate(
        f"选 D={SELECTED_D}\n({SELECTED_D}, {d_value} 次)",
        xy=(SELECTED_D, d_value),
        xytext=(7.1, 25),
        color=INK,
        fontsize=9.5,
        arrowprops={"arrowstyle": "->", "color": SECONDARY, "lw": 1.1},
    )
    ax_false.text(
        15.7,
        FALSE_SLEEP_BUDGET + 1.4,
        f"预算 ≤ {FALSE_SLEEP_BUDGET} 次",
        ha="right",
        color=SECONDARY,
        fontsize=9,
    )
    ax_false.set_title("待睡保持越久，误睡越少", fontsize=12.5, pad=12)
    ax_false.set_xlabel("待睡状态所需静息帧 D")
    ax_false.set_ylabel("误睡次数（63 个短静息样本）")
    ax_false.set_xticks(range(1, 17, 2))
    ax_false.set_ylim(bottom=-2)
    ax_false.yaxis.set_major_locator(MaxNLocator(integer=True))

    # Right: measured maximum wake latency versus W.
    style_axis(ax_wake)
    ax_wake.axvspan(1, SELECTED_W, color=BLUE, alpha=0.08, linewidth=0)
    ax_wake.axhline(
        WAKE_BUDGET,
        color="#d03b3b",
        linestyle=(0, (5, 3)),
        linewidth=1.3,
    )
    ax_wake.plot(
        n,
        wake,
        color=ORANGE,
        linewidth=2.2,
        marker="o",
        markersize=5.5,
        markerfacecolor=ORANGE,
        markeredgecolor=SURFACE,
        markeredgewidth=1.4,
        label="实测/上界 1+W",
    )
    w_value = wake[n.index(SELECTED_W)]
    ax_wake.scatter(
        [SELECTED_W],
        [w_value],
        s=110,
        color=ORANGE,
        edgecolor=INK,
        linewidth=1.1,
        zorder=5,
    )
    ax_wake.annotate(
        f"选 W={SELECTED_W}\n上界 {w_value} 帧 ≈ 50 ms",
        xy=(SELECTED_W, w_value),
        xytext=(4.3, 9.0),
        color=INK,
        fontsize=9.5,
        bbox={"boxstyle": "round,pad=0.25", "fc": SURFACE, "ec": "none", "alpha": 0.88},
        arrowprops={"arrowstyle": "->", "color": SECONDARY, "lw": 1.1},
    )
    ax_wake.text(
        15.7,
        WAKE_BUDGET + 0.35,
        f"延迟预算 {WAKE_BUDGET} 帧（60 Hz）",
        ha="right",
        color="#a02c2c",
        fontsize=9,
    )
    ax_wake.set_title("唤醒确认越久，延迟越高", fontsize=12.5, pad=12)
    ax_wake.set_xlabel("被唤醒中所需保持帧 W")
    ax_wake.set_ylabel("睡着 → 完全醒着的最大延迟（帧）")
    ax_wake.set_xticks(range(1, 17, 2))
    ax_wake.set_ylim(1, 18)
    ax_wake.yaxis.set_major_locator(MaxNLocator(integer=True))

    fig.suptitle(
        "两条曲线分别标定 D 与 W：单计数器 D=W=N 时，误睡达标需 N≥6，"
        "延迟达标需 N≤2",
        x=0.055,
        y=0.985,
        ha="left",
        fontsize=12.5,
        color=INK,
        fontweight="semibold",
    )
    fig.text(
        0.055,
        0.025,
        "结论：四状态把“入睡前防抖 D”和“唤醒 ramp W”拆成两个计数器，取 D=6、W=2。\n"
        "唤醒上界只含 1+W=3 帧；D 不在 ASLEEP → WAKING → AWAKE 路径上。",
        ha="left",
        va="bottom",
        fontsize=9.5,
        color=SECONDARY,
        linespacing=1.45,
    )
    fig.tight_layout(rect=(0, 0.11, 1, 0.93))
    args.png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.png, bbox_inches="tight")
    print(f"wrote {args.png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
