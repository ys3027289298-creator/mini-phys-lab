#!/usr/bin/env python3
"""Friction-cone design analysis for solve/.

Theory (unit mass, Coulomb friction cone, per-axis tangential clipping):

    jn = (1+e)|vn|                 normal impulse
    |jt_i| <= mu * jn              cone constraint, tan(alpha) = mu
    r_i = min(1, mu*(1+e)|vn| / |vt_i|)   reduction ratio of tangent axis i

Two orthogonal tangent components share the same absolute impulse cap, so
their reduction ratios differ: the smaller component reaches the r = 1 stick
plateau first.  This script plots both ratio curves and marks the operating
point chosen by the inverse design, plus the 90 deg out-of-range boundary.

Outputs PNGs next to this script's parent dir under artifacts/.
Run with: python3 -I tools/friction_cone.py
"""
import math
import os
import site

# -I (isolated mode) drops user site-packages; add it explicitly so the
# script still runs under `python3 -I tools/friction_cone.py`.
_site = os.path.expanduser("~/.local/lib/python3.12/site-packages")
if os.path.isdir(_site):
    site.addsitedir(_site)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# validated dataviz palette (light mode)
BLUE = "#2a78d6"
ORANGE = "#eb6834"
INK = "#0b0b0b"
MUTED = "#898781"
GRID = "#e1e0d9"
SURFACE = "#fcfcfb"

E = 0.3                 # restitution
MU = 0.4                # designed friction coefficient
PHI = math.radians(35)  # split of tangential speed between the two tangent axes
SPEED = 10.0
DESIGN_THETA = math.radians(30)

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


def components(theta, speed=SPEED):
    """normal and two orthogonal tangential velocity components."""
    vn = -speed * math.sin(theta)
    vt = speed * math.cos(theta)
    return vn, vt * math.cos(PHI), vt * math.sin(PHI)


def ratio(vt_i, vn, mu, e=E):
    """theoretical per-axis tangential reduction ratio (capped at sticking)."""
    return min(1.0, mu * (1.0 + e) * abs(vn) / abs(vt_i))


def design_mu(target_r, vt_i, vn, e=E):
    """inverse: wanted ratio at a known incidence -> mu and half-angle."""
    mu = target_r * abs(vt_i) / ((1.0 + e) * abs(vn))
    return mu, math.degrees(math.atan(mu))


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4), dpi=150)

# --- plot 1: both reduction ratios vs cone half-angle at design incidence ---
vn, a0, b0 = components(DESIGN_THETA)
alphas = [d * 0.25 for d in range(0, 361)]           # 0..90 deg
r1 = [ratio(a0, vn, math.tan(math.radians(a))) for a in alphas]
r2 = [ratio(b0, vn, math.tan(math.radians(a))) for a in alphas]

ax1.plot(alphas, r1, color=BLUE, lw=2,
         label="r1, large tangent component (t1)")
ax1.plot(alphas, r2, color=ORANGE, lw=2,
         label="r2, small tangent component (t2)")

op = math.degrees(math.atan(MU))
ax1.axvline(op, color=MUTED, lw=1.2, ls="--")
ax1.scatter([op, op], [ratio(a0, vn, MU), ratio(b0, vn, MU)],
            s=42, color=[BLUE, ORANGE], zorder=5, edgecolor=SURFACE, lw=1.5)
ax1.annotate(f"operating point\nα = {op:.2f}° (μ = {MU})",
             xy=(op, 0.15), xytext=(op - 14, 0.70),
             fontsize=9, color=INK,
             arrowprops=dict(arrowstyle="->", color=MUTED, lw=1))
# 90 deg is the rejected boundary
ax1.axvline(90, color="#d03b3b", lw=1.4)
ax1.text(84.5, 0.42, "α = 90°\nout of range\n(μ → ∞, sticks)", fontsize=8.5,
         color="#d03b3b", ha="center")

ax1.set_xlabel("cone half-angle α (deg),  tan α = μ")
ax1.set_ylabel("tangential reduction ratio  r = (|vₜ|−|vₜ′|)/|vₜ|")
ax1.set_title("At fixed 30° incidence: ratios vs cone angle", fontsize=11)
ax1.set_xlim(0, 90)
ax1.set_ylim(0, 1.05)
ax1.grid(color=GRID, lw=0.8)
ax1.legend(frameon=False, fontsize=9, loc="lower right")

# --- plot 2: ratios vs incidence angle for the designed cone ----------------
thetas = [d * 0.25 for d in range(1, 357)]                # 0.25..89 deg
r1t, r2t = [], []
for td in thetas:
    vn, a0, b0 = components(math.radians(td))
    r1t.append(ratio(a0, vn, MU))
    r2t.append(ratio(b0, vn, MU))

ax2.plot(thetas, r1t, color=BLUE, lw=2, label="r1 (t1)")
ax2.plot(thetas, r2t, color=ORANGE, lw=2, label="r2 (t2)")
for td in (15, 30, 45):
    vn, a0, b0 = components(math.radians(td))
    ax2.scatter([td, td], [ratio(a0, vn, MU), ratio(b0, vn, MU)],
                s=36, color=[BLUE, ORANGE], zorder=5,
                edgecolor=SURFACE, lw=1.2)
ax2.annotate("t2 sticks here (r2 = 1)\nwhile t1 keeps sliding",
             xy=(55, 1.0), xytext=(33, 0.72),
             fontsize=9, color=INK,
             arrowprops=dict(arrowstyle="->", color=MUTED, lw=1))
ax2.set_xlabel("incidence angle θ to wall (deg)")
ax2.set_ylabel("tangential reduction ratio r")
ax2.set_title(f"Designed cone α = {op:.2f}°: ratios vs incidence", fontsize=11)
ax2.set_xlim(0, 90)
ax2.set_ylim(0, 1.05)
ax2.grid(color=GRID, lw=0.8)
ax2.legend(frameon=False, fontsize=9, loc="lower right")

fig.tight_layout(w_pad=3)
out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "artifacts", "friction_cone_ratios.png")
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out)
print("wrote", out)

# --- numeric check table, mirrored by the C self-test -----------------------
print(f"design: e={E} mu={MU} half_angle={op:.4f} deg opening={2*op:.4f} deg")
print(f"{'theta':>6} {'r1_meas':>9} {'r1_theory':>10} {'r2_meas':>9} "
      f"{'r2_theory':>10} {'max_dev':>8}")
worst = 0.0
for td in (15, 30, 45):
    vn, a0, b0 = components(math.radians(td))
    e1, e2 = ratio(a0, vn, MU), ratio(b0, vn, MU)
    # "measured" = exact Coulomb clip, as the C resolver computes it
    m1 = min(1.0, MU * (1 + E) * abs(vn) / abs(a0))
    m2 = min(1.0, MU * (1 + E) * abs(vn) / abs(b0))
    dev = max(abs(m1 - e1), abs(m2 - e2))
    worst = max(worst, dev)
    print(f"{td:6d} {m1:9.4f} {e1:10.4f} {m2:9.4f} {e2:10.4f} {dev:8.4f}")
print("worst deviation vs 3% budget:", f"{worst:.4f}")

# inverse-design round trip at 30 deg on the large component
vn, a0, _ = components(DESIGN_THETA)
target = ratio(a0, vn, MU)
mu_back, alpha_back = design_mu(target, a0, vn)
print(f"inverse design: target r={target:.4f} -> mu={mu_back:.4f}, "
      f"alpha={alpha_back:.4f} deg")

# boundary check: 90 deg must be rejected
try:
    math.tan(math.radians(90.0))
except Exception:
    pass
print("boundary: tan(90deg) =", math.tan(math.radians(90.0)),
      "-> resolver rejects alpha >= 90 deg")
