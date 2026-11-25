#!/usr/bin/env python3
"""
Render artifacts/report.html from the experiment JSON.

Pure stdlib: every plotted number is read from artifacts/data, never typed by
hand. Output is a single self-contained HTML file with inline SVG, light/dark
themes, hover tooltips and a full table view (the accessibility twin).
"""
import html
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "artifacts", "data")
OUT = os.path.join(ROOT, "artifacts", "report.html")

TIERS = [("shallow", "浅重叠 (depth 0.15)"),
         ("medium", "中重叠 (depth 0.50)"),
         ("deep", "深重叠 (depth 0.85)")]
COLORS = {"shallow": ("#2a78d6", "#3987e5"),
          "medium": ("#eb6834", "#d95926"),
          "deep": ("#1baf7a", "#199e70")}

W, H = 880, 340
ML, MR, MT, MB = 68, 132, 28, 52
PW, PH = W - ML - MR, H - MT - MB


def load(name):
    with open(os.path.join(DATA, name)) as f:
        return json.load(f)


def sx(x, x0, x1):
    return ML + (x - x0) / (x1 - x0) * PW


def sy(y, y0, y1):
    return MT + PH - (y - y0) / (y1 - y0) * PH


def nice_ticks(lo, hi, n=5):
    step = (hi - lo) / n
    return [lo + i * step for i in range(n + 1)]


def chart_frame(title, subtitle, x_label, y_label, x0, x1, y0, y1,
                xticks, yfmt="{:.0f}".format, body="", log_note=""):
    g = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="{html.escape(title)}">']
    for t in nice_ticks(y0, y1):
        y = sy(t, y0, y1)
        g.append(f'<line x1="{ML}" y1="{y:.1f}" x2="{ML+PW}" y2="{y:.1f}" class="grid"/>')
        g.append(f'<text x="{ML-10}" y="{y+4:.1f}" text-anchor="end" class="tick">{yfmt(t)}</text>')
    for t in xticks:
        x = sx(t, x0, x1)
        g.append(f'<line x1="{x:.1f}" y1="{MT+PH}" x2="{x:.1f}" y2="{MT+PH+6}" class="axis"/>')
        g.append(f'<text x="{x:.1f}" y="{MT+PH+22}" text-anchor="middle" class="tick">{t:g}</text>')
    g.append(f'<line x1="{ML}" y1="{MT+PH}" x2="{ML+PW}" y2="{MT+PH}" class="axis"/>')
    g.append(f'<line x1="{ML}" y1="{MT}" x2="{ML}" y2="{MT+PH}" class="axis"/>')
    g.append(f'<text x="{ML}" y="18" class="title">{html.escape(title)}</text>')
    g.append(f'<text x="{ML}" y="36" class="subtitle">{html.escape(subtitle)}</text>')
    g.append(f'<text x="{ML+PW/2}" y="{H-8}" text-anchor="middle" class="axis-label">{html.escape(x_label)}</text>')
    g.append(f'<text transform="translate(20,{MT+PH/2}) rotate(-90)" text-anchor="middle" class="axis-label">{html.escape(y_label)}</text>')
    if log_note:
        g.append(f'<text x="{ML+PW}" y="{MT-8}" text-anchor="end" class="note">{html.escape(log_note)}</text>')
    g.append(body)
    g.append("</svg>")
    return "".join(g)


def line_series(points, color, tip_prefix):
    pts = " ".join(f"{x:.1f},{y:.1f}" for x, y, _ in points)
    out = [f'<polyline points="{pts}" fill="none" stroke="{color}" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/>']
    for x, y, tip in points:
        out.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4.5" fill="{color}" stroke="var(--surface-1)" stroke-width="2" class="pt" data-tip="{html.escape(tip_prefix + tip)}"/>')
    return "".join(out)


def legend(items, x, y):
    """items: list of (color, label). Lines converge at the right edge, so a
    legend carries identity instead of colliding end labels."""
    out = []
    for k, (color, label) in enumerate(items):
        yy = y + k * 18
        out.append(f'<line x1="{x}" y1="{yy}" x2="{x+22}" y2="{yy}" stroke="{color}" stroke-width="2" stroke-linecap="round"/>')
        out.append(f'<circle cx="{x+11}" cy="{yy}" r="4" fill="{color}" stroke="var(--surface-1)" stroke-width="1.5"/>')
        out.append(f'<text x="{x+28}" y="{yy+4}" class="legend">{html.escape(label)}</text>')
    return "".join(out)


def main():
    runs = load("runs.json")["runs"]
    age = load("age_cdf.json")
    audit = load("frame0_audit.json")
    mem = load("memory.json")
    failures = load("failures.json")
    Hs = sorted({r["H"] for r in runs})
    n_events = sum(r["n_injected"] for r in runs if r["tier"] == "shallow")
    chosen = mem["analytic"]["empirical_choice"]
    h_pick = max(v["min_H_for_1pct"] for v in chosen.values())
    mem_pick = next(c["bytes_24_balls"] for c in mem["curve"] if c["H"] == h_pick)
    mem_1k = next(c["bytes_per_1000_balls"] for c in mem["curve"] if c["H"] == h_pick)

    # ---- Chart 1: recovery rate vs H -------------------------------------
    body = []
    for tier, label in TIERS:
        tr = [r for r in runs if r["tier"] == tier]
        pts = [(sx(r["H"], 0, max(Hs)), sy(r["recovery_rate"], 0, 1),
                f"H={r['H']} 恢复率={r['recovery_rate']*100:.1f}%") for r in tr]
        body.append(line_series(pts, COLORS[tier][0], label + " · "))
    body.append(legend([(COLORS[t][0], l) for t, l in TIERS], ML + PW - 150, MT + 46))
    chart1 = chart_frame(
        "退化恢复率 vs 历史帧数 H",
        "球心重合事件中，由速度或历史帧成功恢复法线的比例（三档初始重叠）",
        "历史帧数 H（帧）", "恢复率", 0, max(Hs), 0, 1, Hs,
        yfmt=lambda v: f"{v*100:.0f}%", body="".join(body))

    # ---- Chart 2: fixed-axis rate in stasis, log-ish linear + 1% line ----
    body = []
    y1 = 1.0
    y_th = sy(0.01, 0, y1)
    body.append(f'<line x1="{ML}" y1="{y_th:.1f}" x2="{ML+PW}" y2="{y_th:.1f}" class="threshold"/>')
    body.append(f'<text x="{ML+PW}" y="{y_th-6:.1f}" text-anchor="end" class="threshold-label">1% 可接受上限</text>')
    for tier, label in TIERS:
        tr = [r for r in runs if r["tier"] == tier]
        pts = []
        for r in tr:
            v = max(r["fixed_in_stasis"], 0.002)   # floor so 0 stays visible
            pts.append((sx(r["H"], 0, max(Hs)), sy(v, 0, y1),
                        f"H={r['H']} 固定轴率={r['fixed_in_stasis']*100:.2f}%"))
        body.append(line_series(pts, COLORS[tier][0], label + " · "))
        cross = next((r["H"] for r in tr if r["fixed_in_stasis"] <= 0.01), None)
        if cross is not None:
            cxp = sx(cross, 0, max(Hs))
            body.append(f'<circle cx="{cxp:.1f}" cy="{y_th:.1f}" r="6" fill="none" stroke="{COLORS[tier][0]}" stroke-width="2"/>')
    body.append(legend([(COLORS[t][0], l) for t, l in TIERS], ML + PW - 150, MT + 46))
    chart2 = chart_frame(
        "静置重合的固定轴兜底率 vs H",
        "仅统计当前速度为零、必须依赖历史窗的事件；零值绘于 0.2% 底线，空心圈为 1% 穿越点",
        "历史帧数 H（帧）", "固定轴兜底率", 0, max(Hs), 0, y1, Hs,
        yfmt=lambda v: f"{v*100:.0f}%", body="".join(body),
        log_note="线性刻度，下限 0.2%")

    # ---- Chart 3: memory vs H --------------------------------------------
    body = []
    pts = [(sx(c["H"], 0, max(Hs)), sy(c["bytes_24_balls"], 0,
            mem["curve"][-1]["bytes_24_balls"] * 1.05),
            f"H={c['H']} 24球={c['bytes_24_balls']}B · 1000球={c['bytes_per_1000_balls']/1024:.0f}KiB")
           for c in mem["curve"]]
    area = f'{ML},{MT+PH} ' + " ".join(f"{x:.1f},{y:.1f}" for x, y, _ in pts) + f' {ML+PW},{MT+PH}'
    body.append(f'<polygon points="{area}" fill="var(--series-1)" opacity="0.10"/>')
    body.append(line_series(pts, "var(--series-1)", "内存 · "))
    xp = sx(h_pick, 0, max(Hs))
    body.append(f'<line x1="{xp:.1f}" y1="{MT}" x2="{xp:.1f}" y2="{MT+PH}" class="threshold"/>')
    body.append(f'<text x="{xp+8:.1f}" y="{MT+14}" class="threshold-label">选定 H={h_pick}</text>')
    chart3 = chart_frame(
        "每对球环形缓冲的内存占用 vs H",
        f"单帧 24 字节（两个 V3）；24 球共 276 对。H={h_pick} 时 {mem_pick/1024:.1f} KiB，1000 球约 {mem_1k/1024/1024:.1f} MiB",
        "历史帧数 H（帧）", "内存（字节）", 0, max(Hs), 0,
        mem["curve"][-1]["bytes_24_balls"] * 1.05, Hs,
        yfmt=lambda v: f"{v/1024:.0f}K", body="".join(body))

    # ---- Chart 4: required-age histogram (pooled over tiers) -------------
    pooled = {}
    for t in ("shallow", "medium", "deep"):
        for row in age[t]["K_histogram"]:
            pooled[row["K"]] = pooled.get(row["K"], 0) + row["count"]
    kmax = max(pooled)
    cmax = max(pooled.values())
    bw = PW / (kmax + 1)
    bars = []
    for k in range(1, kmax + 1):
        cnt = pooled.get(k, 0)
        x = ML + (k - 0.5) * bw
        h = (cnt / cmax) * (PH - 10)
        y = MT + PH - h
        cls = "bar bar-tail" if k >= h_pick else "bar"
        bars.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{max(bw-2,2):.1f}" height="{h:.1f}" rx="4" class="{cls}" data-tip="静置 K={k} 帧 → 需回溯 {k+1} 帧，样本 {cnt}"/>')
    xp = ML + (h_pick - 0.5) * bw
    bars.append(f'<line x1="{xp:.1f}" y1="{MT}" x2="{xp:.1f}" y2="{MT+PH}" class="threshold"/>')
    bars.append(f'<text x="{xp+8:.1f}" y="{MT+14}" class="threshold-label">H={h_pick} 覆盖左侧全部游程</text>')
    xt = [k for k in range(1, kmax + 1, 8)]
    chart4 = chart_frame(
        "静置时长 K 的经验分布（三档汇总）",
        "K 为重合态滞留帧数；恢复所需回溯年龄恰为 K+1，长尾决定 H 的下限",
        "静置时长 K（帧）", "事件数", 1, kmax + 1, 0, cmax, xt,
        yfmt="{:.0f}".format, body="".join(bars))

    # ---- KPI tiles --------------------------------------------------------
    kpis = [
        ("重合事件总数", f"{n_events}", "每档 500 帧 × 0.8 注入率"),
        ("非法线/零向量失效", f"{len(failures)}", "球心重合用例全部解出单位法线"),
        ("选定历史深度", f"H = {h_pick}", f"解析反解 H≥{mem['analytic']['H_for_tail_1pct']}，经验穿越 {h_pick}"),
        (f"H={h_pick} 内存", f"{mem_pick/1024:.1f} KiB", f"24 球；1000 球规模 {mem_1k/1024/1024:.1f} MiB"),
    ]
    kpi_html = "".join(
        f'<div class="tile"><div class="tile-label">{html.escape(l)}</div>'
        f'<div class="tile-value">{html.escape(v)}</div>'
        f'<div class="tile-sub">{html.escape(s)}</div></div>' for l, v, s in kpis)

    # ---- frame-0 audit table ---------------------------------------------
    rows = []
    for a in audit:
        ok = '<span class="ok">✓ 单位法线</span>' if a["ok"] else '<span class="bad">✗ 失效</span>'
        rows.append(
            f"<tr><td>{a['tier']}</td><td>球 {a['ball_a']} / 球 {a['ball_b']}</td>"
            f"<td>{a['kind']}</td><td>{a['source']}</td>"
            f"<td class='num'>({a['normal'][0]:g}, {a['normal'][1]:g})</td>"
            f"<td>{a['history_age']}</td><td>{ok}</td></tr>")
    audit_table = "".join(rows)

    # ---- full run table (accessibility twin) -----------------------------
    full = []
    for r in runs:
        full.append(
            f"<tr><td>{r['tier']}</td><td class='num'>{r['H']}</td>"
            f"<td class='num'>{r['n_injected']}</td>"
            f"<td class='num'>{r['recovery_rate']*100:.2f}%</td>"
            f"<td class='num'>{r['fixed_in_stasis']*100:.2f}%</td>"
            f"<td class='num'>{r['fixed_conditional']*100:.2f}%</td>"
            f"<td class='num'>{r['memory_bytes']}</td></tr>")
    run_table = "".join(full)

    analytic = mem["analytic"]
    doc = f"""<!doctype html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>球心重合退化法线：兜底链与历史帧数定标</title>
<style>
:root {{
  color-scheme: light;
  --surface-1: #fcfcfb; --page: #f9f9f7;
  --text-primary: #0b0b0b; --text-secondary: #52514e; --muted: #898781;
  --grid: #e1e0d9; --axis: #c3c2b7; --border: rgba(11,11,11,0.10);
  --series-1: #2a78d6; --good: #006300; --critical: #d03b3b;
}}
@media (prefers-color-scheme: dark) {{
  :root {{
    color-scheme: dark;
    --surface-1: #1a1a19; --page: #0d0d0d;
    --text-primary: #f0efec; --text-secondary: #c3c2b7; --muted: #898781;
    --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,0.10);
    --series-1: #3987e5; --good: #0ca30c; --critical: #e66767;
  }}
}}
* {{ box-sizing: border-box; }}
body {{ margin:0; background:var(--page); color:var(--text-primary);
  font:14px/1.6 system-ui,-apple-system,"Segoe UI",sans-serif; }}
.wrap {{ max-width:980px; margin:0 auto; padding:32px 24px 64px; }}
h1 {{ font-size:26px; margin:0 0 6px; }}
h2 {{ font-size:18px; margin:36px 0 4px; }}
.lead {{ color:var(--text-secondary); margin:0 0 24px; max-width:70ch; }}
.tiles {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(200px,1fr)); gap:12px; }}
.tile {{ background:var(--surface-1); border:1px solid var(--border); border-radius:12px; padding:16px; }}
.tile-label {{ color:var(--text-secondary); font-size:13px; }}
.tile-value {{ font-size:30px; font-weight:650; margin:2px 0; }}
.tile-sub {{ color:var(--muted); font-size:12px; }}
.card {{ background:var(--surface-1); border:1px solid var(--border); border-radius:12px;
  padding:12px 16px; margin-top:16px; overflow-x:auto; }}
svg {{ width:100%; height:auto; display:block; }}
.grid {{ stroke:var(--grid); stroke-width:1; }}
.axis {{ stroke:var(--axis); stroke-width:1; }}
.tick,.axis-label {{ fill:var(--muted); font-size:11px; font-variant-numeric:tabular-nums; }}
.axis-label {{ font-size:12px; }}
.title {{ fill:var(--text-primary); font-size:14px; font-weight:600; }}
.subtitle,.note,.end-label {{ fill:var(--text-secondary); font-size:11px; }}
.legend {{ fill:var(--text-secondary); font-size:11px; }}
.threshold {{ stroke:var(--critical); stroke-width:1.5; }}
.threshold-label {{ fill:var(--critical); font-size:11px; font-weight:600; }}
.bar {{ fill:var(--series-1); }}
.bar-tail {{ fill:var(--muted); opacity:0.55; }}
.pt {{ cursor:default; }}
table {{ border-collapse:collapse; width:100%; font-size:13px; }}
th,td {{ padding:7px 10px; border-bottom:1px solid var(--grid); text-align:left; }}
th {{ color:var(--text-secondary); font-weight:600; }}
td.num,th.num {{ font-variant-numeric:tabular-nums; }}
.ok {{ color:var(--good); font-weight:600; }}
.bad {{ color:var(--critical); font-weight:600; }}
details {{ margin-top:12px; }}
summary {{ cursor:pointer; color:var(--text-secondary); }}
.callout {{ background:var(--surface-1); border:1px solid var(--border); border-left:4px solid var(--critical);
  border-radius:12px; padding:14px 18px; margin-top:16px; }}
.formula {{ font-family:ui-monospace,Menlo,monospace; background:var(--page);
  border:1px solid var(--border); border-radius:8px; padding:2px 6px; }}
#tip {{ position:fixed; pointer-events:none; background:var(--text-primary); color:var(--surface-1);
  padding:5px 9px; border-radius:7px; font-size:12px; opacity:0; transition:opacity .08s; z-index:9;
  max-width:280px; }}
</style>
</head>
<body>
<div class="wrap">
<h1>球心重合的退化碰撞法线：兜底链与历史帧数定标</h1>
<p class="lead">教科书法线 <span class="formula">n = Δx / |Δx|</span> 在两球球心重合时是 0/0。
本报告先给出极限分析与四级兜底定义，再用三档初始重叠、500 帧 × 24 球的受控实验，
以退化恢复率与内存占用两条曲线反解必须保留的历史帧数。所有法线均由 C 库
<span class="formula">sphere_contact()</span> 经 ctypes 实时计算。</p>

<div class="tiles">{kpi_html}</div>

<h2>1. 极限不确定性与兜底链</h2>
<p class="lead">令 <span class="formula">d(t)=c_B(t)−c_A(t)</span>，法线为 <span class="formula">d/|d|</span>。
球心重合时 d→0：若相对速度 d′ 有限非零，由洛必达法则极限 <span class="formula">d/|d| → d′/|d′|</span> 存在，
但它取决于<strong>到达重合的路径</strong>而非当前位置，几何本身无法确定方向；
若 d′ 也趋于零（静置重合/出生重叠），当前状态不含任何方向信息，只能向历史帧或固定世界轴兜底。</p>
<div class="card">{chart1}</div>
<div class="card">{chart2}</div>

<h2>2. 两条曲线定帧数</h2>
<p class="lead">静置时长 K 服从几何分布（q={analytic['q']}，均值 {analytic['mean_stasis']:.0f} 帧）。
窗口深 H 时取不到证据的概率为 <span class="formula">P(K ≥ H) = (1−q)^H</span>，
令其 ≤ 1% 反解 <span class="formula">H ≥ ln(0.01)/ln(1−q) = {analytic['H_for_tail_1pct']}</span>；
经验曲线在网格点 H={h_pick} 处三档全部穿越 1%。内存随 H 线性增长，每帧每对 24 字节。</p>
<div class="card">{chart4}</div>
<div class="card">{chart3}</div>

<h2>3. 球心重合用例审计（三档各跑一次）</h2>
<p class="lead">每档三个预埋重合对：迎头（速度分支）、静置 2 帧（历史分支）、冷启动空历史（固定轴分支）。
任何崩掉或输出零向量的用例在此点名球编号——实际失效数为 <strong>{len(failures)}</strong>。</p>
<div class="card"><table>
<thead><tr><th>重叠档</th><th>球对</th><th>情形</th><th>法线来源</th><th>法线</th><th>回溯年龄</th><th>结果</th></tr></thead>
<tbody>{audit_table}</tbody></table></div>

<div class="callout"><strong>固定轴兜底的物理代价：只滑不弹。</strong>
固定轴与真实相对运动正交时，法向冲量为零、动能不损失（C 自检
<span class="formula">fixed_axis_no_energy_loss</span> 验证），两球沿轴滑开而非弹开。
因此固定轴只允许作为最后兜底，其在重合事件中的出现率必须低于 1%——H={h_pick} 时三档经验值为
{'; '.join(f"{t} {next(r['fixed_in_stasis'] for r in runs if r['tier']==t and r['H']==h_pick)*100:.2f}%" for t,_ in TIERS)}。</div>

<details><summary>完整实验数据表（33 行，表格视图）</summary>
<div class="card"><table>
<thead><tr><th>重叠档</th><th class="num">H</th><th class="num">事件数</th><th class="num">恢复率</th>
<th class="num">静置固定轴率</th><th class="num">全事件固定轴率</th><th class="num">内存(字节)</th></tr></thead>
<tbody>{run_table}</tbody></table></div></details>
</div>
<div id="tip"></div>
<script>
const tip = document.getElementById('tip');
document.querySelectorAll('.pt, .bar').forEach(el => {{
  el.addEventListener('mousemove', e => {{
    tip.textContent = el.dataset.tip;
    tip.style.opacity = 1;
    tip.style.left = Math.min(e.clientX + 14, window.innerWidth - 300) + 'px';
    tip.style.top = (e.clientY + 14) + 'px';
  }});
  el.addEventListener('mouseleave', () => tip.style.opacity = 0);
}});
</script>
</body></html>"""
    with open(OUT, "w") as f:
        f.write(doc)
    print(f"wrote {OUT} ({len(doc)//1024} KiB)")


if __name__ == "__main__":
    main()
