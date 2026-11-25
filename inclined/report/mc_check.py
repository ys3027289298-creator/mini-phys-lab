#!/usr/bin/env python3
"""Monte-Carlo cross-check of the analytic slip-misjudgement rate.

The C module models the accepted angle as theta_meas ~ N(theta_true, sigma)
and predicts a false-slide rate P(theta_meas > theta_s(mu)) = tail(z).  Here we
sample the same model and report the empirical rate for several m, so the
closed form in inclined.c is checked against actual draws rather than assumed.

Read the binary's own calibration (mu0 = tan 20 deg, sigma = 0.5 deg) so the
two agree by construction on the *model*; the point is that sampling reproduces
1 - Phi((theta_s - 20)/sigma).
"""
import csv
import math
import os
import random
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "out")
BIN = os.path.join(OUT, "inclined")

REQUIRED_DEG = 20.0
SIGMA_DEG = 0.5
N = 200_000


def analytic(mu):
    theta_s = math.degrees(math.atan(mu))
    z = (theta_s - REQUIRED_DEG) / SIGMA_DEG
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def monte_carlo(mu, seed):
    rng = random.Random(seed)
    theta_s = math.degrees(math.atan(mu))
    hit = 0
    for _ in range(N):
        # true angle is exactly 20 deg (the boundary case); measurement noise
        if REQUIRED_DEG + rng.gauss(0.0, SIGMA_DEG) > theta_s:
            hit += 1
    return hit / N


def main():
    rows = []
    with open(os.path.join(OUT, "curves.csv"), newline="") as f:
        for r in csv.DictReader(f):
            m = float(r["m"])
            if abs(m - round(m, 2)) < 1e-9:     # a handful of m values
                rows.append((m, float(r["mu"])))
    print(f"{'m':>6} {'mu':>9} {'analytic':>11} {'MCMC':>11} {'abs_err':>10}")
    worst = 0.0
    for i, (m, mu) in enumerate(rows):
        a = analytic(mu)
        s = monte_carlo(mu, seed=1000 + i)
        err = abs(a - s)
        worst = max(worst, err)
        print(f"{m:6.2f} {mu:9.6f} {a:11.6g} {s:11.6g} {err:10.2e}")
    print(f"worst |analytic - MCMC| = {worst:.3e}  (N={N})")
    # 200k draws put the MC standard error near 1e-3 at p~0.5, so 5e-3 is a
    # generous, sampling-noise-aware tolerance.
    ok = worst < 5e-3
    print(f"mc_matches_closed_form={int(ok)}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
