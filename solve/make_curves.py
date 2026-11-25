#!/usr/bin/env python3
"""Render the measurement figures from artifacts/ CSV + results.json.

Reads  : <repo>/artifacts/{iter_time,convergence,eps_sweep}.csv,results.json
Writes : <repo>/artifacts/fig_*.png
All labels are English so the charts render without a CJK font.
"""
import csv, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ART = os.path.abspath(os.path.join(HERE, "..", "artifacts"))

BLUE, ORANGE, AQUA, RED = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"
INK, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#d8d6d1", "#fcfcfb"

plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.titlecolor": INK, "figure.facecolor": "white", "axes.facecolor": SURF,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.7,
    "axes.spines.top": False, "axes.spines.right": False,
})


def read_csv(name):
    with open(os.path.join(ART, name)) as f:
        return list(csv.DictReader(f))


def fig_time_vs_iters(res):
    rows = read_csv("iter_time.csv")
    K = [int(r["K"]) for r in rows]
    leg = [float(r["legacy_median_us"]) / 1000 for r in rows]
    cac = [float(r["cachedfix_median_us"]) / 1000 for r in rows]
    fl = res["fit_legacy_us"]
    fc = res["fit_cached_us"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.3))

    # left: full range, linear axes -> the relationship is a straight line
    kk = list(range(0, 4100, 50))
    ax1.plot(kk, [(fl["intercept"] + fl["slope_per_iter"] * k) / 1000 for k in kk],
             "--", color=BLUE, lw=1.4, alpha=.8,
             label=f"legacy fit: T={fl['intercept']:.0f}+{fl['slope_per_iter']:.2f}K µs, R²={fl['r2']:.4f}")
    ax1.plot(K, leg, "o", color=BLUE, ms=5, label="legacy (recompute every iter)")
    ax1.plot(kk, [(fc["intercept"] + fc["slope_per_iter"] * k) / 1000 for k in kk],
             "--", color=ORANGE, lw=1.4, alpha=.8,
             label=f"cached fit: slope {fc['slope_per_iter']:.2f} µs/iter, R²={fc['r2']:.4f}")
    ax1.plot(K, cac, "s", color=ORANGE, ms=5, label="cached derived data")
    ax1.axhline(4.0, color=RED, lw=1.8)
    ax1.text(60, 4.15, "budget 4.0 ms", color=RED, fontsize=9, weight="bold")
    ax1.axvline(4096, color=MUTED, ls=":", lw=1.2)
    i = K.index(4096)
    ax1.annotate(f"current K=4096\n{leg[i]:.2f} ms (median)",
                 xy=(4096, leg[i]), xytext=(2350, 9.6),
                 arrowprops=dict(arrowstyle="->", color=INK, lw=1),
                 fontsize=9, weight="bold")
    k4 = fl["K_at_4000us"]
    ax1.plot([k4], [4.0], "D", color=RED, ms=7, mec="white", mew=1.2)
    ax1.annotate(f"4 ms crossing\nK≈{k4:.0f}", xy=(k4, 4.0), xytext=(1150, 2.0),
                 arrowprops=dict(arrowstyle="->", color=RED, lw=1), fontsize=9, color=RED)
    ax1.set_xlabel("fixed iteration count K")
    ax1.set_ylabel("contact-solve time per frame (ms)")
    ax1.set_title("Total time is linear in iterations (30 balls, 106 contacts)")
    ax1.set_ylim(0, 13)
    ax1.legend(fontsize=8.2, loc="upper left", framealpha=.92)

    # right: low-K zoom -> where convergence actually stops
    zK = [k for k in K if k <= 32]
    zleg = leg[:len(zK)]
    zcac = cac[:len(zK)]
    ax2.plot(zK, zleg, "-o", color=BLUE, ms=6, label="legacy")
    ax2.plot(zK, zcac, "-s", color=ORANGE, ms=6, label="cached, fixed K")
    kstar = res["final"]["cold_iters"]
    tstar = res["final"]["median_us"] / 1000
    ax2.axvline(kstar, color=AQUA, ls="--", lw=1.6)
    ax2.plot([kstar], [tstar], "*", color=AQUA, ms=18, mec=INK, mew=.6, zorder=5)
    ax2.annotate(f"ε*=0.01 convergence stop\nK_min={kstar}, measured {tstar*1000:.1f} µs",
                 xy=(kstar, tstar), xytext=(11.5, 0.052),
                 arrowprops=dict(arrowstyle="->", color=AQUA, lw=1.2),
                 fontsize=9, weight="bold")
    ax2.set_xlabel("iteration count K")
    ax2.set_ylabel("ms")
    ax2.set_title("Operating point vs 4096 configured")
    ax2.legend(fontsize=9, loc="upper left")
    ax2.set_ylim(0, max(zleg) * 1.12)

    fig.tight_layout()
    p = os.path.join(ART, "fig_time_vs_iters.png")
    fig.savefig(p, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_convergence():
    rows = read_csv("convergence.csv")
    it = [int(r["iter"]) for r in rows]
    d = [max(float(r["dpn_max"]), 1e-12) for r in rows]
    cut = 40
    it, d = it[:cut], d[:cut]

    fig, ax = plt.subplots(figsize=(8.6, 4.4))
    ax.semilogy(it, d, "-o", color=BLUE, ms=4.5, lw=1.6,
                label=r"$\max_c|\Delta P_n^{(k)}|$ (measured)")
    k3 = next(k for k, dd in zip(it, d) if dd < 0.03)
    k5 = next(k for k, dd in zip(it, d) if dd < 0.01)
    ax.axhline(0.03, color=ORANGE, ls="--", lw=1.5)
    ax.plot([k3], [0.03], "o", color=ORANGE, ms=9, mec="white", mew=1.3, zorder=5)
    ax.annotate(f"loosest pass ε=0.03 → stop at iter {k3}",
                xy=(k3, 0.03), xytext=(9, 0.075),
                fontsize=9.5, color=ORANGE, weight="bold",
                arrowprops=dict(arrowstyle="->", color=ORANGE, lw=1.1))
    ax.axhline(0.01, color=AQUA, ls="--", lw=1.5)
    ax.plot([k5], [0.01], "o", color=AQUA, ms=9, mec="white", mew=1.3, zorder=5)
    ax.annotate(f"chosen ε=0.01 → stop at iter {k5}",
                xy=(k5, 0.01), xytext=(10, 0.0016),
                fontsize=9.5, color=AQUA, weight="bold",
                arrowprops=dict(arrowstyle="->", color=AQUA, lw=1.1))
    ax.text(15.5, 3e-8, r"cold compacted frame: $\Delta P_n(k)\approx0.082\cdot2^{-k}$",
            fontsize=9, color=MUTED)
    ax.set_xlabel("Gauss-Seidel iteration")
    ax.set_ylabel("max normal-impulse change (impulse units, log scale)")
    ax.set_title("Convergence criterion gives the minimum iteration count")
    ax.set_xticks(range(0, 41, 4))
    ax.set_ylim(1e-13, 0.3)
    ax.legend(fontsize=10, loc="lower left")
    fig.tight_layout()
    p = os.path.join(ART, "fig_convergence.png")
    fig.savefig(p, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_threshold_curves():
    rows = read_csv("eps_sweep.csv")
    eps = [float(r["eps"]) for r in rows]
    resmax = [float(r["res_max"]) for r in rows]
    tmed = [float(r["solve_med_us"]) for r in rows]
    tp95 = [float(r["solve_p95_us"]) for r in rows]
    itmed = [int(r["it_med"]) for r in rows]

    fig, (axr, axt, axi) = plt.subplots(
        3, 1, figsize=(8.8, 8.6), sharex=True,
        gridspec_kw={"height_ratios": [2.1, 1.5, 1.1]})

    # residual
    axr.semilogx(eps, resmax, "-o", color=BLUE, ms=6, lw=1.8,
                 label="worst penetration over 2400 frames")
    axr.axhline(1e-3, color=RED, lw=2)
    axr.text(0.9, 1.17e-3,
             "position-residual limit 1e-3  (r = 1, normalized = absolute)",
             color=RED, fontsize=9, weight="bold", ha="left")
    estars = 0.01
    istar = min(range(len(eps)), key=lambda i: abs(eps[i] - estars))
    rstar = resmax[istar]
    axr.axvline(estars, color=AQUA, ls="--", lw=1.6)
    axr.plot([estars], [rstar], "*", color=AQUA, ms=20, mec=INK, mew=.6, zorder=6)
    axr.annotate(f"chosen ε*=0.01\nres_max={rstar:.2e}  (4.5× margin)",
                 xy=(estars, rstar), xycoords="data",
                 xytext=(0.985, 0.96), textcoords="axes fraction",
                 ha="right", va="top",
                 fontsize=10, weight="bold",
                 arrowprops=dict(arrowstyle="->", color=AQUA, lw=1.3))
    ipass = min(range(len(eps)), key=lambda i: abs(eps[i] - 0.04))
    epass, rpass = eps[ipass], resmax[ipass]
    axr.plot([epass], [rpass], "D", color=ORANGE, ms=8, mec="white", mew=1.2)
    axr.annotate("loosest deterministic pass ε≈0.04\n9.87e-4: only 1.3% margin",
                 xy=(epass, rpass), xytext=(0.22, 2.4e-4), fontsize=9, color=ORANGE,
                 arrowprops=dict(arrowstyle="->", color=ORANGE, lw=1.1))
    axr.set_yscale("log")
    axr.set_ylabel("max position residual\n(penetration / r)")
    axr.set_title("Threshold trade-off: position residual and frame time vs impulse threshold",
                  fontsize=14, pad=12)
    axr.legend(fontsize=9, loc="lower left")

    # time
    axt.semilogx(eps, tmed, "-o", color=ORANGE, ms=6, lw=1.8, label="median solve µs")
    axt.semilogx(eps, tp95, ":", color=MUTED, lw=1.4, label="p95 µs")
    axt.axvline(estars, color=AQUA, ls="--", lw=1.4)
    for idx, (e, t, it) in enumerate(zip(eps, tmed, itmed)):
        changed = idx == 0 or it != itmed[idx - 1] or idx == len(eps) - 1
        if changed:
            axt.annotate(f"{it} it", (e, t), textcoords="offset points",
                         xytext=(0, 8), ha="center", fontsize=8.5, color=INK)
    axt.set_ylabel("solve time (µs/frame)")
    axt.set_ylim(0, max(tp95) * 1.22)
    axt.text(0.60, 0.95,
             "4 ms budget = 4000 µs (off scale);\noptimized cold-frame p95 ≤ 0.02 ms",
             transform=axt.transAxes, ha="center", va="top",
             fontsize=8.8, color=MUTED)
    axt.legend(fontsize=8.5, loc="upper left")

    # iterations
    axi.semilogx(eps, itmed, "-s", color=AQUA, ms=6, lw=1.8)
    axi.axvline(estars, color=AQUA, ls="--", lw=1.4)
    axi.set_ylabel("median iters/frame")
    axi.set_xlabel("impulse-change threshold ε (larger = looser  →  smaller = stricter)")
    axi.set_yticks(sorted(set(itmed)))

    for ax in (axr, axt, axi):
        ax.invert_xaxis()
    fig.tight_layout()
    p = os.path.join(ART, "fig_threshold_curves.png")
    fig.savefig(p, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_phase(res):
    b = res["baseline"]
    names = ["contact\nenumeration", "redundant setup\nevery iteration",
             "impulse-apply\ncore"]
    keys = ["enum_us", "redund_us", "core_us"]
    cols = [BLUE, ORANGE, AQUA]
    vals = [b[k] for k in keys]
    shares = [b["enum_share"] * 100, b["redund_share"] * 100, b["core_share"] * 100]

    fig, ax = plt.subplots(figsize=(8.8, 3.9))
    y = range(len(names))
    bars = ax.barh(list(y), vals, color=cols, height=.58,
                   edgecolor="white", linewidth=1.5)
    ax.set_xscale("log")
    ax.set_yticks(list(y))
    ax.set_yticklabels(names)
    ax.invert_yaxis()
    xmax = max(vals) * 4
    ax.set_xlim(1, xmax)
    for yi, v, s in zip(y, vals, shares):
        label = f"{v:,.1f} µs   ({s:.2f}%)" if s >= 0.01 else f"{v:.2f} µs   ({s:.3f}%)"
        ax.text(v * 1.25, yi, label, va="center", fontsize=10, weight="bold")
    ax.set_xlabel(f"median time per {b['outer_median_us']/1000:.2f} ms legacy solve (µs, log scale)")
    ax.set_title("Where the 8 ms goes: timer-measured three-way split (K=4096, 106 contacts)")
    ax.grid(axis="y", visible=False)
    fig.tight_layout()
    p = os.path.join(ART, "fig_phase_breakdown.png")
    fig.savefig(p, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return p


def main():
    with open(os.path.join(ART, "results.json")) as f:
        res = json.load(f)
    for p in [fig_phase(res), fig_time_vs_iters(res), fig_convergence(),
              fig_threshold_curves()]:
        print("wrote", p)


if __name__ == "__main__":
    main()
