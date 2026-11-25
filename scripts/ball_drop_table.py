#!/usr/bin/env python3
"""Ball-drop-to-rest timetable and restitution invariant check.

Model (matches the repo's collide/solve/integrate C modules):
  * point-of-contact height 0, sphere radius r, centre rest height r
  * gravity g, impulse law at impact:  v+ = -e * v-   (0 < e < 1)
  * rigid impacts are instantaneous; flight between impacts is ballistic

Outputs (all under artifacts/ball_drop/):
  bounce_events.csv   analytic impact/peak events
  timetable_10ms.csv  full timetable, one row every 10 ms
  invariant.csv       analytic + sampled peak-height ratio check
  grid_scan.csv       grid-width trade-off data (anomalies vs fill effort)
  grid_tradeoff.html  two-curve chart (anomaly count, fill effort vs grid width)
"""
import argparse
import csv
import math
import os
import time

# ---- physical parameters (SI) -------------------------------------------------
G = 9.8          # gravity, m/s^2 (same constant as integrate.c)
E = 0.70         # coefficient of restitution (solve.c uses v+ = -e v-)
R = 0.05         # sphere radius, m
Y0 = 1.00        # initial centre height at release, m  (z0 = Y0 - R)
DT_MAIN = 0.010  # requested timetable grid, s
DEV_LIMIT = 0.05 # 5 % invariant tolerance

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "artifacts", "ball_drop")


# ---- segment 1/2/3: analytic event construction -------------------------------
def build_events():
    """Return (impacts, peaks, T_inf). impacts[k] = k-th contact (k>=1)."""
    z0 = Y0 - R
    u0 = math.sqrt(2 * G * z0)          # speed at the first impact
    t_f = u0 / G                        # segment 1: free-fall duration

    impacts, peaks = [], []
    t = t_f
    u = u0                              # incoming speed at current impact
    k = 0
    while True:
        k += 1
        # k-th impact: incoming speed u_k = e^(k-1) u0 at time t_k
        impacts.append(dict(k=k, t=t, v_in=-u, v_out=E * u,
                            z=0.0, y=R))
        # k-th bounce: launch speed w_k = e^k u0
        w = E * u
        flight = 2 * w / G              # segment 2 round-trip duration
        apex_t = t + flight / 2
        z_apex = w * w / (2 * G)
        peaks.append(dict(k=k, t=apex_t, t_prev=t, t_next=t + flight,
                          z=z_apex, y=R + z_apex, w=w, flight=flight))
        u = w                           # next incoming speed
        t += flight
        # geometric series is infinite; truncate once physically sub-microsecond
        if flight < 1e-6 and z_apex < 1e-9:
            break
    T_inf = t_f + (2 * E * u0) / (G * (1 - E))   # segment 3: Zeno limit
    return impacts, peaks, T_inf, dict(z0=z0, u0=u0, t_f=t_f)


# ---- analytic state at any instant --------------------------------------------
def analytic_state(t, impacts, peaks, T_inf):
    """Centre height y, vertical velocity vy (up positive), contact marker."""
    if t <= 0.0:
        return Y0, 0.0, 0
    if t >= T_inf:
        return R, 0.0, 2
    if t <= impacts[0]["t"]:                       # segment 1: free fall
        tau = t
        return Y0 - 0.5 * G * tau * tau, -G * tau, 0
    # locate the bounce whose flight window contains t
    for p in peaks:                                # segment 2: k-th round trip
        if p["t_prev"] <= t < p["t_next"]:
            s = t - p["t_prev"]
            w = p["w"]
            z = w * s - 0.5 * G * s * s
            return R + z, w - G * s, 0
    return R, 0.0, 2                               # segment 3: resting


# ---- sampled timetable on a fixed grid ----------------------------------------
def build_grid(dt, impacts, peaks, T_inf):
    """Snap each impact to its bracketing grid row; marker 0 air / 1 hit / 2 rest."""
    n = math.ceil(T_inf / dt) + 1                  # one row beyond T_inf
    hit_rows = {}                                  # row -> analytic impact index k
    for imp in impacts:
        hit_rows.setdefault(round(imp["t"] / dt), imp["k"])
    rows = []
    for i in range(n + 1):
        t = i * dt
        y, vy, marker = analytic_state(t, impacts, peaks, T_inf)
        note = ""
        if i in hit_rows and marker == 0:
            marker = 1                            # grid row samples an impact
            y, vy = R, float("nan")               # velocity flips inside the cell
            note = f"impact #{hit_rows[i]} inside cell: v- -> v+ = -e v-"
        if marker == 2:
            note = "rest (Zeno limit reached)"
        rows.append(dict(i=i, t=t, y=y, vy=vy, contact=marker, note=note))
    return rows, hit_rows


# ---- invariant: analytic peaks and grid-reconstructed peaks -------------------
def check_invariant(dt, rows, impacts, peaks, T_inf):
    # analytic: every pair satisfies H_{k+1}/H_k = e^2 exactly
    analytic = []
    for a, b in zip(peaks, peaks[1:]):
        ratio = b["z"] / a["z"]
        analytic.append(dict(pair=f"{a['k']}-{b['k']}", kind="analytic",
                             k_hi=a["k"],
                             h_hi=a["z"], h_lo=b["z"], ratio=ratio,
                             dev=abs(ratio - E * E) / (E * E), flag=False))

    # sampled: for every grid-resolvable analytic bounce, take the tallest
    # timetable row strictly inside its two bracketing impacts. A bounce is
    # resolvable only when its flight time spans >= 2 grid cells; shorter
    # bounces alias into a cell and are excluded as out-of-resolution, so the
    # Zeno tail cannot manufacture phantom anomalies or duplicate peak labels.
    resolvable = [p for p in peaks if p["flight"] >= 2 * dt]
    impact_row = {}
    for imp in impacts:
        impact_row.setdefault(imp["k"], round(imp["t"] / dt))

    sampled_peaks = []
    for p in resolvable:
        r0 = impact_row[p["k"]]                  # launch impact row
        r1 = round(p["t_next"] / dt)             # next impact row
        interior = [r for r in rows
                    if r0 < r["i"] < r1 and r["contact"] == 0]
        if not interior:
            continue                             # nominally resolvable, but missed
        top = max(interior, key=lambda r: r["y"])
        sampled_peaks.append(dict(k=p["k"], t=top["t"], z=top["y"] - R))

    sampled = []
    for a, b in zip(sampled_peaks, sampled_peaks[1:]):
        ratio = b["z"] / a["z"]
        dev = abs(ratio - E * E) / (E * E)
        sampled.append(dict(pair=f"k{a['k']}->k{b['k']}",
                            kind=f"sampled_{int(dt*1000)}ms", k_hi=a["k"],
                            h_hi=a["z"], h_lo=b["z"], ratio=ratio,
                            dev=dev, flag=dev > DEV_LIMIT))
    n_resolvable = len(resolvable)
    return analytic, sampled, n_resolvable, len(sampled_peaks)


# ---- grid-width trade-off scan ------------------------------------------------
def grid_scan(impacts, peaks, T_inf):
    widths = [0.001, 0.002, 0.005, 0.010, 0.020, 0.025, 0.050, 0.100]
    out = []
    for dt in widths:
        t0 = time.perf_counter()
        rows, _ = build_grid(dt, impacts, peaks, T_inf)
        _, sampled, n_resolvable, n_found = check_invariant(
            dt, rows, impacts, peaks, T_inf)
        wall = time.perf_counter() - t0
        n_bad = sum(1 for s in sampled if s["flag"])
        frontier = min((p["z"] for p in peaks if p["flight"] >= 2 * dt),
                       default=0.0)
        out.append(dict(dt_ms=round(dt * 1000), rows=len(rows),
                        wall_us=round(wall * 1e6, 1),
                        bounces_resolvable=n_resolvable,
                        bounces_seen=n_found,
                        missed_among_resolvable=n_resolvable - n_found,
                        smallest_resolvable_peak_mm=round(frontier * 1000, 2),
                        anomalies_gt5pct=n_bad))
    return out


# ---- self-contained SVG chart (no external deps) ------------------------------
def write_chart(scan, path):
    def x_log(dt_ms):  # log axis 1..100 ms
        return 70 + (math.log10(dt_ms) - 0.0) / 2.0 * 560

    def y_lin(v, vmax, top=60, h=180):
        return top + h - (v / vmax) * h if vmax else top + h

    anom = [s["anomalies_gt5pct"] for s in scan]
    missed = [s["missed_among_resolvable"] for s in scan]
    rows = [s["rows"] for s in scan]
    amax = max(anom + missed) or 1
    rmax = max(rows) or 1

    def poly(vals, vmax, top):
        return " ".join(f"{x_log(s['dt_ms']):.1f},{y_lin(v, vmax, top):.1f}"
                        for s, v in zip(scan, vals))

    def dot(cx, cy, r, color, tip):
        # visible mark plus an >=12px transparent hover target with native tooltip
        return (f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r}" fill="{color}"/>'
                f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="10" fill="transparent">'
                f"<title>{tip}</title></circle>")

    dots_a = "".join(
        (dot(x_log(s["dt_ms"]), y_lin(s["anomalies_gt5pct"], amax), 5, "#eb6834",
             f"格宽 {s['dt_ms']} ms：异常 {s['anomalies_gt5pct']} 对，"
             f"可分辨弹跳 {s['bounces_resolvable']} 个，{s['rows']} 行")
         + f'<text x="{x_log(s["dt_ms"]):.1f}" y="{y_lin(s["anomalies_gt5pct"], amax)-12:.1f}" '
           f'text-anchor="middle" class="lbl">{s["anomalies_gt5pct"]}</text>')
        for s in scan)
    dots_m = "".join(
        dot(x_log(s["dt_ms"]), y_lin(s["missed_among_resolvable"], amax), 4, "#1baf7a",
            f"格宽 {s['dt_ms']} ms：可分辨弹跳漏检 {s['missed_among_resolvable']} 个")
        for s in scan)
    dots_r = "".join(
        dot(x_log(s["dt_ms"]), y_lin(s["rows"], rmax, 300), 5, "#2a78d6",
            f"格宽 {s['dt_ms']} ms：{s['rows']} 行，实测 {s['wall_us']} µs")
        for s in scan)
    xticks = "".join(
        f'<line x1="{x_log(s["dt_ms"]):.1f}" y1="240" x2="{x_log(s["dt_ms"]):.1f}" y2="246" stroke="#94a3b8"/>'
        f'<text x="{x_log(s["dt_ms"]):.1f}" y="262" text-anchor="middle" class="axis">{s["dt_ms"]}</text>'
        f'<line x1="{x_log(s["dt_ms"]):.1f}" y1="480" x2="{x_log(s["dt_ms"]):.1f}" y2="486" stroke="#94a3b8"/>'
        f'<text x="{x_log(s["dt_ms"]):.1f}" y="502" text-anchor="middle" class="axis">{s["dt_ms"]}</text>'
        for s in scan)
    chosen_x = x_log(10)
    html = f"""<!doctype html><html lang="zh"><head><meta charset="utf-8">
<title>格宽取舍曲线</title>
<style>
 body{{font-family:system-ui,sans-serif;margin:24px;color:#0f172a;background:#f8fafc}}
 .card{{background:#fff;border:1px solid #e2e8f0;border-radius:12px;padding:16px 20px;max-width:720px}}
 .axis{{font-size:11px;fill:#475569}} .lbl{{font-size:11px;fill:#475569;font-weight:600}}
 h2{{font-size:15px;margin:4px 0}} p{{font-size:12px;color:#475569}}
</style></head><body><div class="card">
<h2>定格宽：可发现异常次数 与 填表耗时（行数代理）随格宽变化</h2>
<svg viewBox="0 0 700 540" width="700" height="540" role="img" aria-label="grid width tradeoff">
 <line x1="70" y1="60" x2="70" y2="240" stroke="#cbd5e1"/><line x1="70" y1="240" x2="630" y2="240" stroke="#cbd5e1"/>
 <polyline points="{poly(anom, amax, 60)}" fill="none" stroke="#eb6834" stroke-width="2"/>
 {dots_a}
 <polyline points="{poly(missed, amax, 60)}" fill="none" stroke="#1baf7a" stroke-width="2" stroke-dasharray="6 4"/>
 {dots_m}
 <line x1="{chosen_x}" y1="56" x2="{chosen_x}" y2="244" stroke="#94a3b8" stroke-dasharray="4 4" stroke-width="1.5"/>
 <circle cx="86" cy="74" r="5" fill="#eb6834"/><text x="98" y="78" class="axis">偏差&gt;5% 的相邻弹跳对数（异常，实线）</text>
 <circle cx="86" cy="90" r="4" fill="#1baf7a"/><text x="98" y="94" class="axis">可分辨弹跳中的漏检数（虚线）</text>
 <text x="{chosen_x+6}" y="232" class="axis">选定 10 ms</text>
 <text x="350" y="285" text-anchor="middle" class="axis">格宽 (ms, 对数轴)</text>

 <line x1="70" y1="300" x2="70" y2="480" stroke="#cbd5e1"/><line x1="70" y1="480" x2="630" y2="480" stroke="#cbd5e1"/>
 <polyline points="{poly(rows, rmax, 300)}" fill="none" stroke="#2a78d6" stroke-width="2"/>
 {dots_r}
 <line x1="{chosen_x}" y1="296" x2="{chosen_x}" y2="484" stroke="#94a3b8" stroke-dasharray="4 4" stroke-width="1.5"/>
 <circle cx="86" cy="314" r="5" fill="#2a78d6"/><text x="98" y="318" class="axis">时间表行数（填表耗时代理；实测墙钟时间见 grid_scan.csv，悬浮看数值）</text>
 <text x="350" y="525" text-anchor="middle" class="axis">格宽 (ms, 对数轴)</text>
 {xticks}
</svg>
<p>判读：格宽 ≤5 ms 时异常数为 0、漏检为 0，但行数翻倍以上；≥20 ms 后分辨率前沿抬升、可分辨弹跳数下降。
10 ms 是拐点：252 行、无漏检，仅剩 2 对越界且都落在毫米级末段（峰高 3.1 mm / 1.5 mm，已与一格自由落体距离同量级），
更细的 5 ms 才能消除，代价是行数 ×2。按任务指定定格宽 10 ms。</p>
</div></body></html>"""
    with open(path, "w", encoding="utf-8") as f:
        f.write(html)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    impacts, peaks, T_inf, meta = build_events()

    with open(os.path.join(args.out, "bounce_events.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["event", "k", "t_s", "y_m", "z_above_rest_m", "v_m_s"])
        w.writerow(["release", 0, 0.0, Y0, meta["z0"], 0.0])
        for im in impacts:
            w.writerow(["impact", im["k"], f"{im['t']:.6f}", R, 0.0,
                        f"{im['v_in']:.6f} -> {im['v_out']:.6f}"])
        for p in peaks:
            w.writerow(["peak", p["k"], f"{p['t']:.6f}", f"{p['y']:.6f}",
                        f"{p['z']:.6f}", 0.0])
        w.writerow(["rest", "", f"{T_inf:.6f}", R, 0.0, 0.0])

    rows, hit_rows = build_grid(DT_MAIN, impacts, peaks, T_inf)
    with open(os.path.join(args.out, "timetable_10ms.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["row", "t_ms", "centre_height_m", "vy_m_s_up_positive",
                    "contact_0air_1hit_2rest", "note"])
        for r in rows:
            w.writerow([r["i"], round(r["t"] * 1000),
                        f"{r['y']:.4f}",
                        "" if math.isnan(r["vy"]) else f"{r['vy']:.4f}",
                        r["contact"], r["note"]])

    # full Markdown timetable (every 10 ms row, no elision)
    with open(os.path.join(args.out, "timetable_10ms.md"), "w", encoding="utf-8") as f:
        f.write("# 落球 10 ms 时间表（完整，共 %d 行）\n\n" % len(rows))
        f.write("接触标记：0=空中，1=该格内发生接触（冲量瞬时，格内速度由 v⁻ 翻转为 v⁺），2=静止\n\n")
        f.write("| 行 | t (ms) | 球心高度 (m) | 竖直速度 (m/s, 向上+) | 接触 | 备注 |\n")
        f.write("|--:|--:|--:|--:|:--:|:--|\n")
        for r in rows:
            vy = "—（格内翻转）" if math.isnan(r["vy"]) else f"{r['vy']:.4f}"
            f.write(f"| {r['i']} | {round(r['t']*1000)} | {r['y']:.4f} | {vy} | "
                    f"{r['contact']} | {r['note']} |\n")

    analytic, sampled, n_resolvable, n_found = check_invariant(
        DT_MAIN, rows, impacts, peaks, T_inf)
    with open(os.path.join(args.out, "invariant.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pair", "source", "peak_hi_m", "peak_lo_m", "ratio",
                    "target_e2", "deviation_pct", "flag_gt5pct"])
        for d in analytic + sampled:
            w.writerow([d["pair"], d["kind"], f"{d['h_hi']:.6f}",
                        f"{d['h_lo']:.6f}", f"{d['ratio']:.6f}",
                        f"{E*E:.4f}", f"{100*d['dev']:.3f}",
                        "YES" if d["flag"] else "no"])

    scan = grid_scan(impacts, peaks, T_inf)
    with open(os.path.join(args.out, "grid_scan.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(scan[0].keys()))
        w.writeheader()
        w.writerows(scan)
    write_chart(scan, os.path.join(args.out, "grid_tradeoff.html"))

    # console summary
    print(f"z0={meta['z0']:.4f} m  u0={meta['u0']:.4f} m/s  t_f={meta['t_f']:.6f} s")
    print(f"T_inf={T_inf:.6f} s  impacts={len(impacts)}  grid_rows={len(rows)}")
    print(f"resolvable bounces at {int(DT_MAIN*1000)} ms: {n_resolvable}, seen: {n_found}")
    print("grid_scan:")
    for s in scan:
        print("  ", s)
    print("sampled invariant flags:")
    for d in sampled:
        if d["flag"]:
            print("  FLAG", d)
    print("analytic max deviation:", max(d["dev"] for d in analytic))


if __name__ == "__main__":
    main()
