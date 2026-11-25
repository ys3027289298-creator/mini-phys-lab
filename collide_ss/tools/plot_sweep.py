#!/usr/bin/env python3
"""plot_sweep: 读取 sweep.csv，画误判/坏法线两条双曲线（双面板，各自纵轴）。

面板1：全部样本（含大坐标消去噪声）；面板2：普通尺度 L<=10。
eps=0 的崩溃点单独用红色标记；标注交叉点与 1e-4 选定阈值。
"""
import csv
import sys
import os

# -I 隔离模式下用户 site 不在 sys.path；显式补回以便找到 matplotlib，
# 同时保留 -I 对 cwd/脚本目录的隔离（防止下载目录里的冒名模块）。
_usersite = os.path.expanduser("~/.local/lib/python3.12/site-packages")
if os.path.isdir(_usersite) and _usersite not in sys.path:
    sys.path.append(_usersite)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, ORANGE, AQUA, RED = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e1e0d9"
SURFACE = "#fcfcfb"


def load(path):
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            rows.append({k: float(v) for k, v in r.items()})
    return rows


def series(rows, mk, bk):
    """返回 (eps>0 的点, eps=0 的点)。"""
    pos = [r for r in rows if r["eps"] > 0]
    zero = [r for r in rows if r["eps"] == 0]
    return ([r["eps"] for r in pos], [r[mk] for r in pos],
            [r[bk] for r in pos]), (zero[0][mk], zero[0][bk], zero[0]["crashes"])


def panel(ax, rows, mk, bk, title, zero_note):
    (xs, mis, bad), (m0, b0, cr0) = series(rows, mk, bk)
    ax.set_facecolor(SURFACE)
    ax.set_xscale("log")
    ax.set_title(title, fontsize=11, color=INK, loc="left")
    ax.set_xlabel("threshold eps_t  (ball-centre distance)", color=MUTED, fontsize=9)
    ax.set_ylabel("event count over contacts", color=MUTED, fontsize=9)
    ax.grid(True, color=GRID, linewidth=0.8, zorder=0)
    for spine in ax.spines.values():
        spine.set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=8)

    l1, = ax.plot(xs, mis, color=BLUE, lw=2, marker="o", ms=5,
                  label="misjudged (reliable geometry sent to fallback)", zorder=3)
    l2, = ax.plot(xs, bad, color=ORANGE, lw=2, marker="s", ms=5, linestyle="--",
                  label="bad normal (angle >30 deg vs double truth)", zorder=3)

    # eps=0（无守卫）：NaN/坏法线与真崩溃，红色单独标记
    ax.scatter([1e-9], [max(b0, 1)], color=RED, s=90, marker="X", zorder=5,
               label=f"eps=0: {int(cr0)} crashes, {int(b0)} NaN/bad")
    ax.annotate(zero_note, xy=(1e-9, max(b0, 1)),
                xytext=(2.2e-9, max(b0, 1) * 0.45 + 400), color=RED, fontsize=8,
                arrowprops=dict(arrowstyle="->", color=RED, lw=1))

    # 交叉点：两条曲线差距最小的网格点（计数为离散值，严格交点无定义）
    if any(m > 0 for m in mis):
        j = min(range(len(xs)), key=lambda i: abs(mis[i] - bad[i]) if mis[i] > 0 else 1e18)
        ax.axvline(xs[j], color=MUTED, lw=1, linestyle=":")
        ax.annotate(f"crossing grid ~ {xs[j]:.0e}\nmis={int(mis[j])} bad={int(bad[j])}",
                    xy=(xs[j], max(mis[j], bad[j])),
                    xytext=(xs[j] * 1.7, max(mis[j], bad[j]) * 1.5 + 300),
                    fontsize=8, color=MUTED)

    # 选定阈值 1e-4：位于噪声地板(1.2e-7)与交叉点(2e-4)之间，偏保守一侧
    ax.axvline(1e-4, color=AQUA, lw=2, linestyle="-", alpha=0.8, zorder=2)
    ax.annotate("selected eps_t = 1e-4", xy=(1e-4, ax.get_ylim()[1] * 0.92),
                fontsize=8, color=AQUA, rotation=90, va="top", ha="right")

    # 直接标签（CVD 二级编码），不依赖图例颜色
    ax.text(xs[-1], mis[-1], "  misjudged", color=BLUE, fontsize=8, va="center")
    ax.text(xs[-1], max(bad[-1], 1), "  bad normal", color=ORANGE, fontsize=8, va="center")
    ax.legend(handles=[l1, l2], fontsize=7.5, frameon=False, loc="center right")


def main():
    in_path, out_path = sys.argv[1], sys.argv[2]
    rows = load(in_path)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 9.5), facecolor=SURFACE)
    z0 = rows[0]
    note_all = (f"eps=0 no guard\n{int(z0['crashes'])} divide-by-zero crashes\n"
                f"{int(z0['bad_normal'])} NaN/bad normals")
    note_loc = ("eps=0 no guard\nlocal scale: 0 crashes\n(float still resolves d)")
    panel(ax1, rows, "misjudged", "bad_normal",
          "All contacts (mixed local scale + large-coordinate cancellation)", note_all)
    panel(ax2, rows, "misjudged_local", "bad_normal_local",
          "Local scale only (L <= 10): float32 noise floor ~ 1.2e-7, zero errors here",
          note_loc)
    fig.suptitle("Threshold trade-off: misjudged collisions vs corrupted normals\n"
                 "200k sphere-sphere contacts, double-precision truth vs float32 observed",
                 fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_path, dpi=150, facecolor=SURFACE)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
