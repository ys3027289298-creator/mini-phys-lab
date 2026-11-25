#!/usr/bin/env python3
"""Generate a single-file interactive HTML report from broadphase results.csv.

Usage: python3 make_report.py <results.csv> <report.html>
Stdlib only. Run with -I.
"""
import csv
import json
import sys


def load(path):
    tiers = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            r = float(row["radius"])
            key = f"{r:g}"  # stable string key: 1.0 -> "1", matching SUMMARY
            tiers.setdefault(key, []).append({
                "h": float(row["h"]),
                "g": int(row["cells_per_axis"]),
                "cand": int(row["candidates"]),
                "true": int(row["true_pairs"]),
                "miss": int(row["misses"]),
                "fp": int(row["false_pos"]),
                "total": int(row["total_pairs"]),
                "frac": float(row["narrow_frac"]),
                "reduction": float(row["reduction_vs_brute"]),
                "theory": float(row["theory_frac"]),
            })
    for rows in tiers.values():
        rows.sort(key=lambda d: d["h"])
    return tiers


def summarize(tiers):
    out = {}
    for key, rows in tiers.items():
        zero = [d for d in rows if d["miss"] == 0]
        tight = min(zero, key=lambda d: d["h"])
        first_bad = min((d for d in rows if d["miss"] > 0),
                        key=lambda d: d["h"], default=None)
        out[key] = {
            "tight_h": tight["h"],
            "tight_cand": tight["cand"],
            "tight_red": tight["reduction"],
            "first_bad_h": first_bad["h"] if first_bad else None,
            "first_bad_miss": first_bad["miss"] if first_bad else None,
            "true_total": rows[-1]["true"],
        }
    return out


HTML = """<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>宽相盒尺寸取舍 — mini-phys-lab broadphase</title>
<style>
.viz-root {
  color-scheme: light;
  --surface-1: #fcfcfb; --page: #f9f9f7;
  --text-primary: #0b0b0b; --text-secondary: #52514e; --muted: #898781;
  --grid: #e1e0d9; --axis: #c3c2b7; --border: rgba(11,11,11,0.10);
  --s1: #2a78d6; --s2: #eb6834; --s3: #1baf7a;
  --theory: #898781; --threshold: #d03b3b;
}
@media (prefers-color-scheme: dark) {
  :root:where(:not([data-theme="light"])) .viz-root {
    color-scheme: dark;
    --surface-1: #1a1a19; --page: #0d0d0d;
    --text-primary: #f0efec; --text-secondary: #c3c2b7; --muted: #898781;
    --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,0.10);
    --s1: #3987e5; --s2: #d95926; --s3: #199e70;
    --theory: #9a988f; --threshold: #e66767;
  }
}
:root[data-theme="dark"] .viz-root {
  color-scheme: dark;
  --surface-1: #1a1a19; --page: #0d0d0d;
  --text-primary: #f0efec; --text-secondary: #c3c2b7; --muted: #898781;
  --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,0.10);
  --s1: #3987e5; --s2: #d95926; --s3: #199e70;
  --theory: #9a988f; --threshold: #e66767;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--page); color: var(--text-primary);
  font-family: system-ui, -apple-system, "Segoe UI", sans-serif; line-height: 1.6; }
.viz-root { max-width: 1080px; margin: 0 auto; padding: 32px 28px 64px; }
h1 { font-size: 26px; margin: 0 0 4px; }
h2 { font-size: 19px; margin: 36px 0 10px; }
h3 { font-size: 15px; margin: 20px 0 6px; color: var(--text-secondary); }
p, li { font-size: 14px; color: var(--text-secondary); }
.sub { color: var(--muted); font-size: 13px; margin-top: 0; }
.math { background: var(--surface-1); border: 1px solid var(--border);
  border-radius: 8px; padding: 14px 18px; margin: 12px 0; font-size: 14px;
  color: var(--text-primary); overflow-x: auto; }
.math .eq { display: block; margin: 6px 0; white-space: nowrap; }
sup, sub { font-size: 0.72em; }
.chart-card { background: var(--surface-1); border: 1px solid var(--border);
  border-radius: 12px; padding: 18px 18px 8px; margin: 16px 0; }
.chart-title { font-size: 15px; font-weight: 600; margin: 0 0 2px; }
.chart-note { font-size: 12px; color: var(--muted); margin: 0 0 10px; }
svg { width: 100%; height: auto; display: block; }
.grid-line { stroke: var(--grid); stroke-width: 1; }
.axis-line { stroke: var(--axis); stroke-width: 1; }
.tick-text { fill: var(--muted); font-size: 11px; }
.axis-label { fill: var(--text-secondary); font-size: 12px; }
.series { fill: none; stroke-width: 2; }
.dot { stroke: var(--surface-1); stroke-width: 1.5; }
.end-label { font-size: 12px; font-weight: 600; dominant-baseline: middle; }
.vline { stroke-width: 1.5; stroke-dasharray: 5 4; }
.hline { stroke-width: 1.5; stroke-dasharray: 6 4; }
.ref-label { font-size: 11px; fill: var(--muted); }
.crosshair { stroke: var(--muted); stroke-width: 1; stroke-dasharray: 3 3; }
.tooltip { position: fixed; pointer-events: none; background: var(--surface-1);
  border: 1px solid var(--border); border-radius: 8px; padding: 8px 10px;
  font-size: 12px; box-shadow: 0 4px 16px rgba(0,0,0,0.18); display: none;
  z-index: 10; min-width: 150px; }
.tooltip .tt-title { font-weight: 600; margin-bottom: 4px; color: var(--text-primary); }
.tooltip .tt-row { display: flex; justify-content: space-between; gap: 12px;
  color: var(--text-secondary); font-variant-numeric: tabular-nums; }
.tt-swatch { display: inline-block; width: 9px; height: 9px; border-radius: 2px;
  margin-right: 5px; }
table { border-collapse: collapse; width: 100%; font-size: 12.5px;
  font-variant-numeric: tabular-nums; background: var(--surface-1);
  border: 1px solid var(--border); border-radius: 10px; overflow: hidden; }
th, td { padding: 6px 10px; text-align: right; border-bottom: 1px solid var(--grid); }
th:first-child, td:first-child { text-align: left; }
th { color: var(--text-secondary); font-weight: 600; background: var(--page);
  position: sticky; top: 0; }
td { color: var(--text-secondary); }
tr.miss td { color: var(--threshold); font-weight: 600; }
.legend { display: flex; flex-wrap: wrap; gap: 16px; margin: 8px 0 0;
  font-size: 12.5px; color: var(--text-secondary); }
.legend .sw { display: inline-block; width: 14px; height: 3px; border-radius: 2px;
  vertical-align: middle; margin-right: 6px; }
.kpi-row { display: grid; grid-template-columns: repeat(auto-fit, minmax(230px, 1fr));
  gap: 12px; margin: 14px 0; }
.kpi { background: var(--surface-1); border: 1px solid var(--border);
  border-radius: 10px; padding: 12px 14px; }
.kpi .k-label { font-size: 12px; color: var(--muted); }
.kpi .k-value { font-size: 22px; font-weight: 700; margin: 2px 0; }
.kpi .k-sub { font-size: 12px; color: var(--text-secondary); }
.theme-toggle { float: right; font-size: 12px; }
.theme-toggle button { font: inherit; padding: 4px 10px; border-radius: 6px;
  border: 1px solid var(--border); background: var(--surface-1);
  color: var(--text-secondary); cursor: pointer; }
.table-wrap { max-height: 420px; overflow: auto; border-radius: 10px; }
code { background: var(--page); border: 1px solid var(--border); border-radius: 4px;
  padding: 1px 5px; font-size: 12.5px; }
</style>
</head>
<body>
<div class="viz-root" id="root">
<div class="theme-toggle">
  <button onclick="document.documentElement.setAttribute('data-theme','light')">浅色</button>
  <button onclick="document.documentElement.setAttribute('data-theme','dark')">深色</button>
</div>
<h1>宽相包围盒剔除：盒尺寸 h 的取舍</h1>
<p class="sub">场景 N=400 个球，L=100×100 平面均匀分布；三档等半径 r ∈ {0.5, 1.0, 2.0}，每档固定随机种子。全量两两检测 = N(N−1)/2 = 79,800 次。</p>

<h2>1. 为什么盒相交后仍必须跑窄相</h2>
<p>两个球的包围盒（以球心为中心、边长 2r 的正方形；网格实现里等价为“所在单元 + 3×3 邻域”）相交只说明球心坐标差在各轴上都不超过阈值，并不约束欧氏距离。对角相邻时，盒在角点接触而球心相距 √2·h 量级，可远大于 r<sub>i</sub>+r<sub>j</sub>。因此宽相输出的是<strong>候选对</strong>，每个候选对仍需一次精确距离判据：</p>
<div class="math">
  <span class="eq">窄相命中 ⟺ d<sub>ij</sub><sup>2</sup> = (x<sub>i</sub>−x<sub>j</sub>)<sup>2</sup> + (y<sub>i</sub>−y<sub>j</sub>)<sup>2</sup> ≤ (r<sub>i</sub>+r<sub>j</sub>)<sup>2</sup></span>
</div>
<p>实验里 <code>false_pos</code> 一列就是“盒留存但窄相拒绝”的假阳性：在 r=0.5、h=8 时 4,339 个候选里 4,318 个是假阳性，仅 21 个真碰撞。</p>

<h2>2. 盒剔除后仍需窄相检测的比例（由 r/h 与 h/L 决定）</h2>
<p>均匀分布、每单元平均 λ = N(h/L)<sup>2</sup> 个球。每个球扫描 3×3 共 9 个单元，候选对数量与窄相留存比例（网格数 g=L/h ≫ 1 时）：</p>
<div class="math">
  <span class="eq">C<sub>bp</sub> ≈ 9N·λ/2 = (9/2)·N<sup>2</sup>·(h/L)<sup>2</sup></span>
  <span class="eq">f<sub>narrow</sub> = C<sub>bp</sub> / C<sub>all</sub> ≈ 9·(h/L)<sup>2</sup>，其中 C<sub>all</sub> = N(N−1)/2</span>
</div>
<p>f<sub>narrow</sub> 只依赖 <strong>h/L</strong>，与球半径 r 无关：假阳性的尺度完全由盒尺寸决定。r 只通过零漏检约束 h ≥ 2r 出现在另一侧（见第 4 节）。g ≤ 3 的粗网格下边界单元邻域不足 9 格，实测值低于该式（角格只扫 4 格），图 B 中可以看到实测曲线与理论虚线在细网格段重合、粗网格段分离。</p>

<h2>3. 用比例反解盒尺寸</h2>
<p>给定可接受的窄相预算 f<sup>*</sup>（占全量两两检测的比例上限），由上式反解：</p>
<div class="math">
  <span class="eq">h ≤ L·√(f<sup>*</sup>/9) = (L/3)·√f<sup>*</sup></span>
  <span class="eq">例：f<sup>*</sup> = 1% → h ≤ 0.0333L；f<sup>*</sup> = 50% → h ≤ 0.236L</span>
</div>

<h2>4. 零漏检约束与最终取值区间</h2>
<p>碰撞球心间距 ≤ r<sub>i</sub>+r<sub>j</sub>。要保证两球所在单元索引之差不超过 1，必须 h ≥ r<sub>i</sub>+r<sub>j</sub>；等半径下即 <strong>h ≥ 2r</strong>。与预算约束联立，盒尺寸的可行窗口为：</p>
<div class="math">
  <span class="eq">2·r<sub>max</sub> ≤ h ≤ (L/3)·√f<sup>*</sup></span>
</div>
<p>窗口为空（2r<sub>max</sub> &gt; L√f<sup>*</sup>/3）说明单一均匀网格无法同时满足零漏检与预算，需要分层网格或扫掠-排序等方案。取窗口下沿 h=2r 即“最紧零漏检盒”，窄相调用最少。</p>

<div id="kpis" class="kpi-row"></div>

<h2>5. 两条曲线定盒尺寸</h2>
<div class="chart-card">
  <p class="chart-title">图 A — 漏检数 vs 盒尺寸 h</p>
  <p class="chart-note">横轴 log 尺度。竖虚线标出各档的 h=2r；越过该线向左漏检迅速上升，向右（含该线）漏检恒为 0。</p>
  <div id="chart-a"></div>
  <div class="legend">
    <span><span class="sw" style="background:var(--s1)"></span>r = 0.5（h=2r=1）</span>
    <span><span class="sw" style="background:var(--s2)"></span>r = 1.0（h=2r=2）</span>
    <span><span class="sw" style="background:var(--s3)"></span>r = 2.0（h=2r=4）</span>
  </div>
</div>
<div class="chart-card">
  <p class="chart-title">图 B — 窄相检测占全量两两检测的比例 vs 盒尺寸 h</p>
  <p class="chart-note">实线为实测 f_narrow，灰色虚线为理论 9(h/L)²（g≫1 成立）；红色点线为 f*=50% 失效阈值。</p>
  <div id="chart-b"></div>
  <div class="legend">
    <span><span class="sw" style="background:var(--s1)"></span>r = 0.5</span>
    <span><span class="sw" style="background:var(--s2)"></span>r = 1.0</span>
    <span><span class="sw" style="background:var(--s3)"></span>r = 2.0</span>
    <span><span class="sw" style="background:var(--theory)"></span>理论 9(h/L)²</span>
    <span><span class="sw" style="background:var(--threshold)"></span>f* = 0.5 阈值</span>
  </div>
</div>

<h2>6. 退化情形：r ≪ h 时宽相失效的判据</h2>
<p>球半径远小于盒尺寸时，单元内球数 λ = N(h/L)² 很大，3×3 邻域几乎覆盖全场景，宽相退化为全量两两检测。由于 f<sub>narrow</sub> ≈ 9(h/L)² 与 r 无关，<strong>失效判据也只由 h/L 决定</strong>：</p>
<div class="math">
  <span class="eq">失效（f<sub>narrow</sub> ≥ f<sup>*</sup>）⟺ h/L ≥ √(f<sup>*</sup>/9)</span>
  <span class="eq">完全失效（f<sub>narrow</sub> → 1）⟺ h/L ≳ 1/3，等价于 g = L/h ≲ 3</span>
  <span class="eq">占用率形式：λ = N(h/L)² ≥ Nf<sup>*</sup>/9（本场景 N=400、f*=0.5 → λ ≳ 22 球/单元）</span>
</div>
<p>实测：h=32（g=3）时 f<sub>narrow</sub>=0.59，h=45.3（g=2）时所有单元互为邻域，f<sub>narrow</sub>=1.00——宽相完全不起作用，与判据一致。此时真碰撞仅数十对，假阳性占候选的 99.5% 以上。</p>

<h2>7. 原始数据（可滚动表格）</h2>
<div class="table-wrap">
<table id="data-table">
<thead><tr>
  <th>r</th><th>h</th><th>g=L/h</th><th>候选数</th><th>真碰撞</th><th>漏检</th>
  <th>假阳性</th><th>窄相占比</th><th>相对全量下降</th><th>理论占比</th>
</tr></thead>
<tbody></tbody>
</table>
</div>
<p class="sub">数据由 broadphase/src/broadphase.c 输出（make sweep），全量基线 79,800 次/帧。</p>
</div>

<div class="tooltip" id="tooltip"></div>
<script>
const TIERS = __DATA__;
const SUMMARY = __SUMMARY__;
const COLORS = ["var(--s1)", "var(--s2)", "var(--s3)"];
const RADII = Object.keys(TIERS).map(Number).sort((a,b)=>a-b);

/* ---------- KPI cards ---------- */
document.getElementById("kpis").innerHTML = RADII.map((r, i) => {
  const s = SUMMARY[r];
  return `<div class="kpi">
    <div class="k-label"><span class="tt-swatch" style="background:${COLORS[i]}"></span>r = ${r}：最紧零漏检点 h = 2r = ${s.tight_h}</div>
    <div class="k-value">${(s.tight_red*100).toFixed(2)}%</div>
    <div class="k-sub">窄相 ${s.tight_cand} 次 / 79,800（下降 ${(s.tight_red*100).toFixed(2)}%），漏检 0；真碰撞 ${s.true_total} 对</div>
  </div>`;
}).join("");

/* ---------- table ---------- */
const tbody = document.querySelector("#data-table tbody");
RADII.forEach(r => TIERS[r].forEach(d => {
  const tr = document.createElement("tr");
  if (d.miss > 0) tr.className = "miss";
  tr.innerHTML = `<td>${r}</td><td>${d.h.toFixed(3)}</td><td>${d.g}</td>
    <td>${d.cand}</td><td>${d.true}</td><td>${d.miss}</td><td>${d.fp}</td>
    <td>${(d.frac*100).toFixed(3)}%</td><td>${(d.reduction*100).toFixed(3)}%</td>
    <td>${(d.theory*100).toFixed(3)}%</td>`;
  tbody.appendChild(tr);
}));

/* ---------- chart helpers ---------- */
const W = 980, H = 380, M = {t: 18, r: 96, b: 48, l: 62};
const iw = W - M.l - M.r, ih = H - M.t - M.b;
const xMin = 0.25, xMax = 128;
const lx = h => M.l + (Math.log10(h) - Math.log10(xMin)) /
                  (Math.log10(xMax) - Math.log10(xMin)) * iw;
const ly = (v, yMax) => M.t + ih - (v / yMax) * ih;
const css = name => getComputedStyle(document.getElementById("root"))
                       .getPropertyValue(name).trim();

function frame(svgId, yMax, yTicks, yFmt) {
  const NS = "http://www.w3.org/2000/svg";
  const host = document.getElementById(svgId);
  host.innerHTML = "";
  const svg = document.createElementNS(NS, "svg");
  svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
  svg.setAttribute("width", W);
  svg.setAttribute("height", H);
  host.appendChild(svg);
  const el = (tag, attrs) => {
    const e = document.createElementNS(NS, tag);
    for (const k in attrs) e.setAttribute(k, attrs[k]);
    svg.appendChild(e); return e;
  };
  yTicks.forEach(v => {
    el("line", {class: "grid-line", x1: M.l, x2: M.l+iw, y1: ly(v,yMax), y2: ly(v,yMax)});
    const t = el("text", {class: "tick-text", x: M.l-8, y: ly(v,yMax)+4,
                          "text-anchor": "end"});
    t.textContent = yFmt(v);
  });
  const hTicks = [0.25,0.5,1,2,4,8,16,32,64,128];
  hTicks.forEach(h => {
    const t = el("text", {class: "tick-text", x: lx(h), y: M.t+ih+18,
                          "text-anchor": "middle"});
    t.textContent = h;
    el("line", {class: "axis-line", x1: lx(h), x2: lx(h),
                y1: M.t+ih, y2: M.t+ih+5});
  });
  el("line", {class: "axis-line", x1: M.l, x2: M.l+iw, y1: M.t+ih, y2: M.t+ih});
  el("line", {class: "axis-line", x1: M.l, x2: M.l, y1: M.t, y2: M.t+ih});
  const xl = el("text", {class: "axis-label", x: M.l+iw/2, y: H-6,
                         "text-anchor": "middle"});
  xl.textContent = "盒尺寸 h（对数轴）";
  return {svg, el};
}

function addSeries(ctx, rows, yMax, access, color, label) {
  const pts = rows.map(d => `${lx(d.h)},${ly(access(d), yMax)}`).join(" ");
  ctx.el("polyline", {class: "series", stroke: color, points: pts});
  rows.forEach(d => ctx.el("circle", {class: "dot", cx: lx(d.h),
      cy: ly(access(d), yMax), r: 3.5, fill: color}));
}

function attachCrosshair(svgId, rows, yMax, tooltipRows) {
  const svg = document.getElementById(svgId);
  const tip = document.getElementById("tooltip");
  const ch = document.createElementNS("http://www.w3.org/2000/svg","line");
  ch.setAttribute("class","crosshair"); ch.setAttribute("y1", M.t);
  ch.setAttribute("y2", M.t+ih); ch.setAttribute("visibility","hidden");
  svg.appendChild(ch);
  svg.addEventListener("mousemove", ev => {
    const rect = svg.getBoundingClientRect();
    const px = (ev.clientX - rect.left) / rect.width * W;
    if (px < M.l || px > M.l+iw) { tip.style.display = "none";
      ch.setAttribute("visibility","hidden"); return; }
    const logH = Math.log10(xMin) + (px-M.l)/iw *
                 (Math.log10(xMax)-Math.log10(xMin));
    const h = Math.pow(10, logH);
    ch.setAttribute("x1", px); ch.setAttribute("x2", px);
    ch.setAttribute("visibility","visible");
    let html = `<div class="tt-title">h ≈ ${h.toFixed(2)}</div>`;
    tooltipRows.forEach((tr, i) => {
      const d = rows[RADII[i]].reduce((best, cur) =>
          Math.abs(Math.log(cur.h)-Math.log(h)) <
          Math.abs(Math.log(best.h)-Math.log(h)) ? cur : best);
      html += `<div class="tt-row"><span><span class="tt-swatch"
        style="background:${COLORS[i]}"></span>${tr.label}</span>
        <span>${tr.fmt(d)}</span></div>`;
    });
    tip.innerHTML = html;
    tip.style.display = "block";
    tip.style.left = Math.min(ev.clientX+14, window.innerWidth-190)+"px";
    tip.style.top = (ev.clientY+14)+"px";
  });
  svg.addEventListener("mouseleave", () => {
    tip.style.display = "none"; ch.setAttribute("visibility","hidden");
  });
}

/* ---------- chart A: misses ---------- */
{
  const yMax = 400;
  const ctx = frame("chart-a", yMax, [0,100,200,300,400], v => v);
  RADII.forEach((r, i) => {
    addSeries(ctx, TIERS[r], yMax, d => d.miss, COLORS[i], `r=${r}`);
    const x = lx(2*r);
    ctx.el("line", {class: "vline", x1: x, x2: x, y1: M.t, y2: M.t+ih,
                    stroke: COLORS[i]});
    const t = ctx.el("text", {class: "ref-label", x: x+4, y: M.t+12});
    t.textContent = `h=2r=${2*r}`;
  });
  const yl = ctx.el("text", {class:"axis-label", x: 16, y: M.t+ih/2,
      "text-anchor":"middle", transform:`rotate(-90 16 ${M.t+ih/2})`});
  yl.textContent = "漏检碰撞对数";
  attachCrosshair("chart-a", TIERS, yMax,
      RADII.map(r => ({label:`r=${r} 漏检`, fmt: d => d.miss})));
}

/* ---------- chart B: narrow fraction ---------- */
{
  const yMax = 1.0;
  const ctx = frame("chart-b", yMax, [0,0.25,0.5,0.75,1],
                    v => (v*100).toFixed(0)+"%");
  const theoryPts = [];
  for (let i=0;i<=64;i++) {
    const h = xMin*Math.pow(xMax/xMin, i/64);
    const f = Math.min(1, 9*(h/100)*(h/100));
    theoryPts.push(`${lx(h)},${ly(f,yMax)}`);
  }
  ctx.el("polyline", {class:"series", stroke: css("--theory"),
      "stroke-dasharray":"7 5", points: theoryPts.join(" ")});
  ctx.el("line", {class:"hline", x1: M.l, x2: M.l+iw,
      y1: ly(0.5,yMax), y2: ly(0.5,yMax), stroke: css("--threshold")});
  const ft = ctx.el("text", {class:"ref-label", x: M.l+8,
      y: ly(0.5,yMax)-5, fill: css("--threshold")});
  ft.textContent = "f* = 0.5 失效阈值";
  RADII.forEach((r, i) =>
      addSeries(ctx, TIERS[r], yMax, d => d.frac, COLORS[i], `r=${r}`));
  const tt = ctx.el("text", {class:"ref-label", x: lx(64), y: ly(9*(64/100)**2,yMax)-8});
  tt.textContent = "理论 9(h/L)²";
  const yl = ctx.el("text", {class:"axis-label", x: 16, y: M.t+ih/2,
      "text-anchor":"middle", transform:`rotate(-90 16 ${M.t+ih/2})`});
  yl.textContent = "窄相调用 / 全量两两检测";
  attachCrosshair("chart-b", TIERS, yMax,
      RADII.map(r => ({label:`r=${r} 占比`,
        fmt: d => (d.frac*100).toFixed(2)+"%"})));
}
</script>
</body>
</html>
"""


def main():
    tiers = load(sys.argv[1])
    summary = summarize(tiers)
    html = (HTML.replace("__DATA__", json.dumps(tiers))
                .replace("__SUMMARY__", json.dumps(summary)))
    with open(sys.argv[2], "w") as f:
        f.write(html)
    for key, s in summary.items():
        print(f"r={key}: tightest zero-miss h={s['tight_h']} "
              f"candidates={s['tight_cand']} reduction={s['tight_red']:.6f} "
              f"first_bad_h={s['first_bad_h']} first_bad_miss={s['first_bad_miss']}")


if __name__ == "__main__":
    main()
