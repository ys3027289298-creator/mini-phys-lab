#!/usr/bin/env python3
"""Substep count trade-off for a ball crossing a thin floor.

The two competing bounds are:

    anti-tunnelling:      s_max(N) < h
    anti-saturation:      s_min(N) > eps

For a frame duration T split into N equal substeps, a monotone fall ending at
speed V has largest semi-implicit Euler displacement s_max ~= V*T/N.  Starting
from (or resting at) zero vertical speed, the first gravity-only substep has
s_min = g*(T/N)^2.  These give an integer interval for N.

The script sweeps N for three frame durations, counts worst-phase discrete
misses against a sphere-vs-slab band h + 2r, and counts position updates that
are snapped to zero by a displacement truncation/dead-zone eps.  It writes a
CSV and a two-metric curve figure under integrate/artifacts/.

Run from the integrate/ directory with:
    python3 -I tools/substep_study.py
or:
    make plot
"""

import csv
import math
import os
import site
import sys

# `python3 -I` ignores user site-packages; add the normal install location so a
# user-installed matplotlib remains usable in isolated mode.
_user_site = os.path.expanduser("~/.local/lib/python3.12/site-packages")
if os.path.isdir(_user_site):
    site.addsitedir(_user_site)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Validated light-mode dataviz palette.
BLUE = "#2a78d6"
ORANGE = "#eb6834"
GOOD = "#0ca30c"
CRITICAL = "#d03b3b"
INK = "#0b0b0b"
MUTED = "#898781"
GRID = "#e1e0d9"
SURFACE = "#fcfcfb"

G = 9.81                 # m/s^2
Y0 = 5.0                 # initial centre height above the floor top, m
V_IMPACT = math.sqrt(2.0 * G * Y0)
H = 0.10                 # floor slab thickness, m
R = 0.02                 # ball radius for the feasible main case, m
EPS = 1.0e-4             # displacement below which movement is truncated to zero, m
RUN_SECONDS = 1.20       # includes impact and a short resting tail
PHASES = 401             # phase offsets in the worst-case miss sweep
ZERO_LIMIT = 0.05
FRAME_DTS = (1.0 / 60.0, 1.0 / 30.0, 0.10)
PROBE_N = (8, 16, 32)    # deliberately over-split probes, one per frame tier
SWEEP_MAX_N = 40

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(os.path.dirname(HERE), "artifacts")


def n_safety_min(dt, speed=V_IMPACT, h=H):
    """Smallest N for speed*dt/N < h under strict inequality."""
    return math.floor(speed * dt / h) + 1


def n_live_max_speed(dt, speed=V_IMPACT, eps=EPS):
    """Largest N for the representative step speed*dt/N > eps."""
    return math.ceil(speed * dt / eps) - 1


def n_live_max_gravity(dt, g=G, eps=EPS):
    """Largest N for the zero-start step g*(dt/N)^2 > eps.

    Solves eps*N^2 < g*dt^2; strict inequality gives ceil(root)-1.
    """
    root = dt * math.sqrt(g / eps)
    return math.ceil(root) - 1


def phase_misses(dt, n, radius=R, h=H, phases=PHASES):
    """Count phase offsets whose two endpoint tests both miss the slab.

    A sphere of radius r against slab [-h, 0] overlaps when its centre y is in
    the expanded band [-h-r, r], width h+2r.  A downward step of length s can
    skip that band iff s > h+2r; the sweep counts how many uniform impact
    phases actually have both endpoints outside it.
    """
    s = V_IMPACT * dt / n
    lower = -h - radius
    upper = radius
    misses = 0
    for j in range(phases):
        phase = (j + 0.5) / phases * s
        y_before = upper + phase
        y_after = y_before - s
        hit_before = lower <= y_before <= upper
        hit_after = lower <= y_after <= upper
        if not hit_before and not hit_after and y_after < lower:
            misses += 1
    return misses


def fall_run(dt, n, eps=EPS):
    """Run one deterministic fall; return tunnel and zero-displacement counts."""
    tau = dt / n
    frames = round(RUN_SECONDS / dt)
    y = Y0
    vy = 0.0
    tunnels = 0
    zero_steps = 0
    contact_steps = 0
    max_step = 0.0

    for _ in range(frames):
        for _ in range(n):
            y_before = y
            vy -= G * tau
            raw_dy = vy * tau

            # Explicit displacement truncation/dead-zone.  Velocity is still
            # integrated, which is why an over-split resting ball accumulates v
            # and then emits a quantized move instead of moving every substep.
            if abs(raw_dy) < eps:
                dy = 0.0
                zero_steps += 1
            else:
                dy = raw_dy
            y += dy
            max_step = max(max_step, abs(dy))

            crossed = y_before > R and y < -H - R
            if crossed:
                tunnels += 1

            hit = -H - R <= y <= R
            if hit:
                contact_steps += 1
                y = R
                if vy < 0.0:
                    vy = 0.0

    total_steps = frames * n
    return {
        "dt": dt,
        "n": n,
        "frames": frames,
        "total_steps": total_steps,
        "zero_steps": zero_steps,
        "zero_ratio": zero_steps / total_steps,
        "tunnels": tunnels,
        "contact_steps": contact_steps,
        "max_step": max_step,
        "final_y": y,
    }


def first_zero_penetration_n(dt):
    for n in range(1, SWEEP_MAX_N + 1):
        if phase_misses(dt, n) == 0:
            return n
    return None


def write_csv(rows):
    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, "substep_sweep.csv")
    fields = [
        "frame_dt", "fps", "n_substeps", "substep_dt",
        "phase_misses", "zero_steps", "total_steps", "zero_ratio",
        "n_safety_min", "n_gravity_max",
    ]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    return path


def make_plot(table):
    plt.rcParams.update({
        "font.family": ["Ubuntu Sans", "DejaVu Sans", "sans-serif"],
        "text.color": INK,
        "axes.edgecolor": "#c3c2b7",
        "axes.labelcolor": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "axes.facecolor": SURFACE,
        "figure.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
    })

    fig, axes = plt.subplots(3, 2, figsize=(11.2, 9.2), dpi=160, sharex="col")
    ns = list(range(1, SWEEP_MAX_N + 1))

    for row, dt in enumerate(FRAME_DTS):
        data = table[dt]
        misses = [data[n]["misses"] for n in ns]
        ratios = [100.0 * data[n]["run"]["zero_ratio"] for n in ns]
        n_min = n_safety_min(dt)
        n_max = n_live_max_gravity(dt)
        n_geom = first_zero_penetration_n(dt)

        ax_p = axes[row][0]
        ax_z = axes[row][1]

        for ax in (ax_p, ax_z):
            span_label = "feasible integer window" if row == 0 and ax is ax_p else None
            ax.axvspan(n_min, min(n_max, SWEEP_MAX_N), color=GOOD, alpha=0.10,
                       label=span_label)
            ax.grid(color=GRID, lw=0.8)
            ax.set_xlim(1, SWEEP_MAX_N)

        ax_p.plot(ns, misses, color=BLUE, lw=2, marker="o", ms=3.2,
                  label="missed phases" if row == 0 else None)
        ax_p.axvline(n_min, color=GOOD, lw=1.6, ls="--",
                     label=f"h-bound N≥{n_min}" if row == 0 else None)
        if n_geom is not None and n_geom != n_min:
            ax_p.axvline(n_geom, color=MUTED, lw=1.2, ls=":",
                         label=f"sphere-band zero at N={n_geom}" if row == 0 else None)
        ax_p.set_ylabel("missed phases / 401")
        ax_p.set_title(f"T = {1000*dt:.2f} ms ({1/dt:.0f} fps): penetration curve",
                       fontsize=10.5)
        if row == 2:
            ax_p.set_xlabel("substeps per frame N")

        ax_z.plot(ns, ratios, color=ORANGE, lw=2, marker="o", ms=3.2,
                  label="zero-displacement steps" if row == 0 else None)
        ax_z.axhline(100.0 * ZERO_LIMIT, color=CRITICAL, lw=1.6, ls="--",
                     label="5% limit" if row == 0 else None)
        ax_z.axvline(n_max, color=INK, lw=1.2, ls=":",
                     label=f"gravity bound N≤{n_max}" if row == 0 else None)
        ax_z.set_ylim(-0.5, 22)
        ax_z.set_ylabel("zero steps (%)")
        ax_z.set_title(f"T = {1000*dt:.2f} ms ({1/dt:.0f} fps): truncation curve",
                       fontsize=10.5)
        if row == 2:
            ax_z.set_xlabel("substeps per frame N")

        probe_n = PROBE_N[row]
        probe_ratio = 100.0 * data[probe_n]["run"]["zero_ratio"]
        if probe_ratio > 5.0:
            ax_z.scatter([probe_n], [probe_ratio], s=58, color=CRITICAL,
                         zorder=5, edgecolor=SURFACE, linewidth=1.4)
            ax_z.annotate(f"N={probe_n}: {probe_ratio:.1f}%\nuse N≤{n_max}",
                          xy=(probe_n, probe_ratio),
                          xytext=(max(2, probe_n - 15), 18.0),
                          fontsize=8.6, color=INK,
                          arrowprops=dict(arrowstyle="->", color=MUTED, lw=1))

    handles0, labels0 = axes[0][0].get_legend_handles_labels()
    handles1, labels1 = axes[0][1].get_legend_handles_labels()
    fig.legend(handles0 + handles1, labels0 + labels1,
               loc="upper center", ncol=4, frameon=False, fontsize=8.8,
               bbox_to_anchor=(0.5, 0.995))
    fig.tight_layout(rect=(0, 0, 1, 0.965), h_pad=2.2, w_pad=2.0)
    path = os.path.join(OUT_DIR, "substep_curves.png")
    fig.savefig(path)
    plt.close(fig)
    return path


def main():
    print(f"V_impact={V_IMPACT:.4f} m/s  h={H:.3f} m  r={R:.3f} m  eps={EPS:.1e} m")
    print(f"strict all-step timestep interval: "
          f"{math.sqrt(EPS/G)*1000:.3f} ms < tau < {H/V_IMPACT*1000:.3f} ms")
    threshold_eps = G * H * H / (V_IMPACT * V_IMPACT)
    print(f"empty-interval threshold eps >= g*h^2/V^2 = {threshold_eps:.4e} m")
    print()

    table = {}
    csv_rows = []
    final_runs = []
    probe_runs = []

    for tier, dt in enumerate(FRAME_DTS):
        table[dt] = {}
        n_min = n_safety_min(dt)
        n_speed_max = n_live_max_speed(dt)
        n_g_max = n_live_max_gravity(dt)
        n_geom = first_zero_penetration_n(dt)

        for n in range(1, SWEEP_MAX_N + 1):
            misses = phase_misses(dt, n)
            run = fall_run(dt, n)
            table[dt][n] = {"misses": misses, "run": run}
            csv_rows.append({
                "frame_dt": f"{dt:.10f}",
                "fps": round(1.0 / dt),
                "n_substeps": n,
                "substep_dt": f"{dt/n:.8f}",
                "phase_misses": misses,
                "zero_steps": run["zero_steps"],
                "total_steps": run["total_steps"],
                "zero_ratio": f"{run['zero_ratio']:.8f}",
                "n_safety_min": n_min,
                "n_gravity_max": n_g_max,
            })

        chosen = n_min
        final = table[dt][chosen]["run"]
        final_phase = table[dt][chosen]["misses"]
        final_runs.append((dt, chosen, final_phase, final))

        probe_n = PROBE_N[tier]
        probe = table[dt][probe_n]["run"]
        probe_phase = table[dt][probe_n]["misses"]
        probe_runs.append((dt, probe_n, probe_phase, probe))

        n5 = max(n for n in range(1, SWEEP_MAX_N + 1)
                 if table[dt][n]["run"]["zero_ratio"] <= ZERO_LIMIT)
        print(f"T={1000*dt:6.2f} ms: N_safe>={n_min:2d}, "
              f"sphere-band zero N={n_geom}, speed bound N<={n_speed_max:4d}, "
              f"gravity/5% bound N<={n_g_max:2d} (scan N5={n5})")

    csv_path = write_csv(csv_rows)
    plot_path = make_plot(table)

    print("\nchosen minimum-call confirmation runs:")
    print(f"{'T(ms)':>8} {'N':>3} {'tau(ms)':>8} {'phase_miss':>10} "
          f"{'tunnels':>8} {'zero/total':>12} {'zero%':>7} {'max_dy(m)':>10}")
    all_ok = True
    for dt, n, phase_miss, run in final_runs:
        ok = phase_miss == 0 and run["tunnels"] == 0 and run["zero_ratio"] <= ZERO_LIMIT
        all_ok &= ok
        print(f"{1000*dt:8.2f} {n:3d} {1000*dt/n:8.3f} {phase_miss:10d} "
              f"{run['tunnels']:8d} {run['zero_steps']:5d}/{run['total_steps']:<5d} "
              f"{100*run['zero_ratio']:6.2f}% {run['max_step']:10.5f} "
              f"{'PASS' if ok else 'FAIL'}")

    print("\nover-split probes (>5% must be corrected downward):")
    for dt, n, phase_miss, run in probe_runs:
        n_max = n_live_max_gravity(dt)
        bad = run["zero_ratio"] > ZERO_LIMIT
        print(f"T={1000*dt:6.2f} ms N={n:2d}: tunnels={run['tunnels']}, "
              f"zero={100*run['zero_ratio']:5.2f}% -> "
              f"{'OVER 5%; choose N<=' + str(n_max) if bad else 'within 5%'}")

    print(f"\nwrote {csv_path}")
    print(f"wrote {plot_path}")
    if not all_ok:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
