#!/usr/bin/env python3
"""Smoke test for measure_video.py: render a synthetic squashing ball clip."""
import os
import subprocess
import sys
import tempfile

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
W, H, FPS = 320, 240, 60
R_PX, BALL_RGB = 28, (227, 73, 72)
N_FRAMES, DROP = 60, 120.0
G_PX = 2.0 * DROP / (0.5 * N_FRAMES / FPS) ** 2   # land at mid-clip


def main():
    work = tempfile.mkdtemp(prefix="smoke_measure_")
    frames = os.path.join(work, "frames")
    os.makedirs(frames)
    floor = H - 12
    t_land = 0.5 * N_FRAMES / FPS
    for i in range(1, N_FRAMES + 1):
        t = i / FPS
        if t <= t_land:
            cy = floor - R_PX - 0.5 * G_PX * (t_land - t) ** 2
            squash = 0.0
        else:  # contact phase: cosine-shaped squash pulse
            tau = (t - t_land) / (0.12)
            squash = 0.15 * max(0.0, np.cos(np.pi * min(tau, 1.0)))
            cy = floor - R_PX * (1.0 - squash)
        sy = 1.0 - squash
        sx = 1.0 / np.sqrt(sy)
        img = Image.new("RGB", (W, H), (252, 252, 251))
        yy, xx = np.mgrid[0:H, 0:W]
        inside = ((xx - W / 2) / (R_PX * sx)) ** 2 + \
                 ((yy - cy) / (R_PX * sy)) ** 2 <= 1.0
        arr = np.asarray(img).copy()
        arr[inside] = BALL_RGB
        Image.fromarray(arr).save(os.path.join(frames, f"f{i:04d}.png"))

    out_csv = os.path.join(work, "measurements.csv")
    r = subprocess.run(
        [sys.executable, os.path.join(HERE, "measure_video.py"), frames,
         "--rgb", *map(str, BALL_RGB), "--tol", "0.25", "--out", out_csv],
        capture_output=True, text=True)
    print(r.stdout)
    if r.returncode != 0:
        print(r.stderr, file=sys.stderr)
        sys.exit("measure_video.py failed on synthetic clip")
    peak = max(float(line.split("epsilon=")[1])
               for line in r.stdout.splitlines() if "epsilon=" in line)
    if abs(peak - 0.15) > 0.02:
        sys.exit(f"expected epsilon ~0.15, measured {peak}")
    print(f"smoke_ok measured_peak={peak:.4f} (truth 0.1500)")


if __name__ == "__main__":
    main()
