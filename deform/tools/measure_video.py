#!/usr/bin/env python3
"""
measure_video.py — assisted frame-by-frame squash measurement.

Given a directory of frames extracted from the reference clip, e.g.

  ffmpeg -i ref.mp4 -vsync 0 frames/f%04d.png

it thresholds the ball by colour, measures its vertical pixel span D on
every frame, and reports the maximum squash ratio per impact:

  epsilon = (D_rest - D_min) / D_rest

D_rest is the median span over the airborne (round, uncompressed) frames.
Pass the ball colour as --rgb R G B and a tolerance --tol 0..1 (fraction
of the RGB channel range); tune it on one frame first with --show-range.

Outputs:
  * per-frame CSV (frame, centre_x, centre_y, width, height, squash)
  * one row per detected impact: impact frame, D_min, epsilon

The ball colour is unknown until the real clip arrives — this script is
shipped unrun.  Smoke test: tools/smoke_measure.sh builds a synthetic
bouncing-ball clip and feeds it through.
"""
import argparse
import csv
import glob
import os
import sys

import numpy as np
from PIL import Image


def ball_mask(arr, rgb, tol):
    target = np.asarray(rgb, dtype=np.float32)
    span = np.abs(arr.astype(np.float32) - target).max(axis=2)
    return span <= tol * 255.0


def measure_frame(path, rgb, tol):
    arr = np.asarray(Image.open(path).convert("RGB"))
    mask = ball_mask(arr, rgb, tol)
    ys, xs = np.nonzero(mask)
    if ys.size < 20:
        return None
    return {
        "cx": 0.5 * (xs.min() + xs.max()),
        "cy": 0.5 * (ys.min() + ys.max()),
        "w": xs.max() - xs.min() + 1,
        "h": ys.max() - ys.min() + 1,
        "pixels": ys.size,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("frames_dir", help="directory of numbered PNG frames")
    ap.add_argument("--rgb", nargs=3, type=int, required=True, metavar=("R", "G", "B"))
    ap.add_argument("--tol", type=float, default=0.25)
    ap.add_argument("--out", default="measurements.csv")
    args = ap.parse_args()

    frames = sorted(glob.glob(os.path.join(args.frames_dir, "f*.png")))
    if not frames:
        sys.exit(f"no frames under {args.frames_dir} (expected f*.png)")

    rows = []
    for idx, path in enumerate(frames, start=1):
        m = measure_frame(path, args.rgb, args.tol)
        if m is None:
            rows.append({"frame": idx, "cx": "", "cy": "", "w": "",
                         "h": "", "squash": ""})
            continue
        rows.append({"frame": idx, "cx": round(m["cx"], 2),
                     "cy": round(m["cy"], 2), "w": m["w"], "h": m["h"],
                     "squash": ""})

    heights = np.asarray([r["h"] for r in rows if r["h"] != ""], dtype=float)
    if heights.size < 10:
        sys.exit("ball detected on too few frames; check --rgb/--tol")
    d_rest = float(np.median(heights))
    for r in rows:
        if r["h"] != "":
            r["squash"] = round(max(0.0, (d_rest - r["h"]) / d_rest), 4)

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["frame", "cx", "cy", "w", "h", "squash"])
        w.writeheader()
        w.writerows(rows)

    # impacts = local minima of vertical span separated by airborne frames
    hs = np.asarray([r["squash"] if r["squash"] != "" else 0.0 for r in rows])
    print(f"D_rest={d_rest:.1f}px  frames={len(frames)}  wrote {args.out}")
    print("impacts (frame of max squash):")
    i = 1
    n = len(hs)
    while i < n - 1:
        if hs[i] > hs[i - 1] and hs[i] >= hs[i + 1] and hs[i] > 0.02:
            j = i
            while j + 1 < n and hs[j + 1] == hs[i]:
                j += 1
            print(f"  impact frame {rows[i]['frame']:4d}: D_min="
                  f"{rows[i]['h']:4d}px  epsilon={hs[i]:.4f}")
            i = j + 1
        else:
            i += 1


if __name__ == "__main__":
    main()
