#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
energy_ledger.py — 落球-弹起能量账 / 粒度实验

四账项（全部用带符号的“能量变化贡献”记账，损失为负、注入为正）：
    N  法向恢复项   N = -1/2 m [(vn-)^2 - (vn+)^2] = -1/2 m (vn-)^2 (1-e^2)
    T  切向摩擦项   T = -1/2 m [(vt-)^2 - (vt+)^2]
                       vt+ = max(|vt-| - mu |vn-|, 0)
    C  数值修正项   C = + mg * dy_corr + 1/2 m (|v_surg+|^2 - |v_surg-|^2)
                       （位置投影注入势能；速度手术带走/注入动能，带符号）
    R  残差（反解） R = dE_obs - N - T - C

于是恒等式 N + T + C + R == dE_obs 严格成立（浮点误差量级）。
“闭合”有两道闸：
    delta_id      = |N+T+C+R-dE| / E0        恒等闭合（构造上即成立）
    unexplained   = |R| / E0                 实质性闭合：残差占初始能量比例
四档落高要求两者都 < 1%（delta_id 实际在 1e-12 量级）。

脚本直接运行会：
  1) 跑四档静止落高，写 four_height_ledger.csv/.json；
  2) 跑 8 个注入缺陷，做五档记账粒度实验，写 granularity_*.csv；
  3) 对每个不闭合/故障账户用 worst_discrepancy() 挑出差额最大的一项；
  4) 出两张图：四档能量瀑布、粒度双曲线（含开销面板）。
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from dataclasses import asdict, dataclass
from typing import Dict, List, Optional

GRAV = 9.81

TERM_CN = {
    "N": "法向恢复",
    "T": "切向摩擦",
    "C": "数值修正",
    "R": "残差",
}

# 数据视觉规范（dataviz 校验通过的 adjacent 顺序）
TERM_COLOR = {
    "N": "#2a78d6",   # blue
    "T": "#eb6834",   # orange
    "C": "#1baf7a",   # aqua
    "R": "#6250d6",   # violet
}
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
SURFACE = "#fcfcfb"


# --------------------------------------------------------------------------- #
# 仿真与记账核心
# --------------------------------------------------------------------------- #
@dataclass
class Params:
    m: float = 1.0
    e: float = 0.5
    mu: float = 0.2
    dt: float = 1e-3
    g: float = GRAV
    beta: float = 1.0     # 位置修正比例（Baumgarte percent），仓库现状 = 1
    slop: float = 0.0     # 穿透容差


def energy(z: float, vx: float, vy: float, p: Params,
           include_tangent: bool = True) -> float:
    """机械能。z 是相对接触位置的球心高度（z=0 恰好接触）。"""
    tx = vx * vx if include_tangent else 0.0
    return p.m * p.g * z + 0.5 * p.m * (tx + vy * vy)


def simulate(
    H: float = 4.0,
    vx0: float = 0.0,
    p: Optional[Params] = None,
    # —— 记账/实现缺陷开关（粒度实验用）——
    speed_N: bool = False,        # D1: 把 |vn|-|vn'|（速度差）当能量入账
    no_T: bool = False,           # D2/D3: 切向摩擦项记 0
    trunc_E: bool = False,        # D3: 总能量漏算 vx^2（只记竖直口径）
    no_C: bool = False,           # D4: 位置修正完全不入账
    e_book: Optional[float] = None,   # D6: 账本用错恢复系数
    double_C: bool = False,       # D7: 修正实际执行两次，账本只记一次
    reverse_friction: bool = False,   # D8: 摩擦不封项，可反向
    max_seconds: float = 30.0,
) -> Dict:
    """落到第一次弹起后的最高点为止，逐帧留下分段计量记录。"""
    if p is None:
        p = Params()
    if e_book is None:
        e_book = p.e

    z, vx, vy = float(H), float(vx0), 0.0
    Efull0 = energy(z, vx, vy, p, True)
    Eobs0 = energy(z, vx, vy, p, not trunc_E)

    frames: List[Dict] = []
    contact_seen = False
    contact_frame = -1
    n_steps = int(max_seconds / p.dt) + 16

    for k in range(n_steps):
        rec: Dict = {"k": k, "t": k * p.dt}
        E_start_obs = energy(z, vx, vy, p, not trunc_E)
        rec["z_pre"] = z

        # ---- 阶段 1：半隐式欧拉积分（与 integrate.c 同一更新顺序） ----
        vy_before_int = vy
        vy -= p.g * p.dt
        z += vy * p.dt
        E_after_int_obs = energy(z, vx, vy, p, not trunc_E)
        dE_int = E_after_int_obs - E_start_obs

        # ---- 阶段 2：碰撞冲量（与 solve.c 同一法向/切向模型） ----
        contact = z < 0.0
        N_book = T_book = 0.0
        N_true = T_true = 0.0
        N_formula = T_formula = 0.0
        vn_pre = vt_pre = vn_post = vt_post = 0.0
        red = 0.0
        reversed_tangent = False
        dE_solve = 0.0
        Jn = Jt = 0.0

        if contact:
            vn_pre, vt_pre = vy, vx
            if vn_pre < 0.0:
                vn_post = -p.e * vn_pre
                u = abs(vt_pre)
                if reverse_friction:
                    red = p.mu * abs(vn_pre)          # 故意不 min(..., u)
                else:
                    red = min(p.mu * abs(vn_pre), u)  # 正确封项
                denom = u if u > 1e-12 else 1.0
                scale = (u - red) / denom
                if u <= 1e-12:
                    scale = 0.0
                vx_new = scale * vt_pre
                reversed_tangent = (vt_pre * vx_new < -1e-9)

                E_pre_solve_obs = energy(z, vx, vy, p, not trunc_E)
                vy = vn_post
                vx = vx_new
                E_post_solve_obs = energy(z, vx, vy, p, not trunc_E)
                dE_solve = E_post_solve_obs - E_pre_solve_obs

                # 独立的状态差口径（真值）
                N_true = -0.5 * p.m * (vn_pre * vn_pre - vn_post * vn_post)
                T_true = -0.5 * p.m * (vt_pre * vt_pre - vx * vx)
                # 系数表式口径
                N_formula = -0.5 * p.m * vn_pre * vn_pre * (1.0 - e_book * e_book)
                T_formula = -(p.m * u * red - 0.5 * p.m * red * red)

                if speed_N:
                    # 仓库 solve.c 现有 out_normal_loss 的量纲：速度差 (m/s)
                    N_book = -p.m * (abs(vn_pre) - abs(vn_post))
                else:
                    N_book = N_formula
                T_book = 0.0 if no_T else T_formula

                Jn = p.m * (vn_post - vn_pre)
                Jt = p.m * (vx - vt_pre)

        # ---- 阶段 3：穿透位置修正 + 速度手术（与 integrate.c 对应） ----
        # 注意：integrate.c 里 y+=pen 与 “if vy<0 vy=0” 都在 if(pen>0) 块内，
        # 速度手术只允许发生在穿透帧，不能每帧清零重力速度。
        pen = -z if z < 0.0 else 0.0
        dy_logged = p.beta * max(pen - p.slop, 0.0) if pen > 0.0 else 0.0
        dy_actual = dy_logged + (dy_logged if double_C else 0.0)
        z += dy_actual

        vy_pre_clamp = vy
        clamped = False
        if pen > 0.0 and vy < 0.0:
            vy = 0.0
            clamped = True
        dK_clamp = 0.5 * p.m * (vy * vy - vy_pre_clamp * vy_pre_clamp)

        E_end_obs = energy(z, vx, vy, p, not trunc_E)
        dE_corr = E_end_obs - E_after_int_obs - dE_solve if contact or dE_solve else \
            E_end_obs - E_after_int_obs
        # 更稳妥：直接重算修正阶段前后
        E_before_corr_obs = E_after_int_obs + dE_solve
        dE_corr = E_end_obs - E_before_corr_obs

        C_true = p.m * p.g * dy_actual + dK_clamp
        C_formula_logged = p.m * p.g * dy_logged + dK_clamp
        if no_C:
            C_book = 0.0
        else:
            C_book = C_formula_logged

        dE_total = E_end_obs - E_start_obs
        R_book = dE_total - N_book - T_book - C_book

        rec.update({
            "contact": bool(contact and vn_pre < 0.0),
            "pen": pen, "dy_logged": dy_logged, "dy_actual": dy_actual,
            "clamped": clamped,
            "z_post": z, "vx": vx, "vy": vy,
            "vn_pre": vn_pre, "vt_pre": vt_pre,
            "vn_post": vn_post, "vt_post": vx if contact else vt_pre,
            "Jn": Jn, "Jt": Jt,
            "reversed_tangent": reversed_tangent,
            "dE_int": dE_int, "dE_solve": dE_solve, "dE_corr": dE_corr,
            "dK_clamp": dK_clamp, "dE": dE_total,
            "N_book": N_book, "T_book": T_book, "C_book": C_book, "R": R_book,
            "N_true": N_true, "T_true": T_true, "C_true": C_true,
            "N_formula": N_formula, "T_formula": T_formula,
            "E_end_obs": E_end_obs,
        })
        frames.append(rec)

        if rec["contact"] and not contact_seen:
            contact_seen = True
            contact_frame = k

        # 第一次弹起后的最高点：vy 由正转非正的那一帧结束
        if contact_seen and k > contact_frame and vy <= 0.0:
            break
        if z < -1.0 * max(H, 1.0):  # 兜底：穿地失控
            break

    # ---- 汇总（直接求和，与 Eend-Estart 互为冗余校验） ----
    def s(key: str) -> float:
        return float(sum(f[key] for f in frames))

    N, T, C, R = s("N_book"), s("T_book"), s("C_book"), s("R")
    dE = s("dE")
    Efull_end = energy(z, vx, vy, p, True)

    out = {
        "params": asdict(p),
        "flags": {"speed_N": speed_N, "no_T": no_T, "trunc_E": trunc_E,
                  "no_C": no_C, "double_C": double_C,
                  "reverse_friction": reverse_friction,
                  "e_book": p.e if e_book is None else e_book},
        "H": H, "vx0": vx0,
        "n_frames": len(frames),
        "contact_frame": contact_frame,
        "z_apex": z, "vx_apex": vx, "vy_apex": vy,
        "E0": Efull0,
        "E_apex_full": Efull_end,
        "Eobs0": Eobs0,
        "Eobs_end": frames[-1]["E_end_obs"] if frames else Eobs0,
        "dE": dE, "N": N, "T": T, "C": C, "R": R,
        "N_true": s("N_true"), "T_true": s("T_true"), "C_true": s("C_true"),
        "frames": frames,
    }
    # 独立的全状态口径残差（不受 trunc_E 影响）：
    # 用真实机械能变化 + 真值三项反解
    out["R_fullstate"] = (Efull_end - Efull0) - out["N_true"] - out["T_true"] - out["C_true"]
    return out


# --------------------------------------------------------------------------- #
# 闭合检查与“差额最大项”挑选
# --------------------------------------------------------------------------- #
def closure_metrics(run: Dict) -> Dict:
    E0 = run["E0"]
    ident = abs(run["N"] + run["T"] + run["C"] + run["R"] - run["dE"])
    unexpl = abs(run["R"])
    return {
        "identity_err_pctE0": 100.0 * ident / E0,
        "residual_pctE0": 100.0 * unexpl / E0,
        "identity_err_pct_dE": 100.0 * ident / max(abs(run["dE"]), 1e-30),
        "pass_1pct": (100.0 * ident / E0 < 1.0) and (100.0 * unexpl / E0 < 1.0),
    }


def term_discrepancies(run: Dict) -> Dict[str, float]:
    """每个账项与其独立直接计量值的差额（J）。残差项取 |R|。"""
    return {
        "N": abs(run["N"] - run["N_true"]),
        "T": abs(run["T"] - run["T_true"]),
        "C": abs(run["C"] - run["C_true"]),
        "R": abs(run["R"]),
    }


def worst_discrepancy(run: Dict, exact_frac: float = 1e-7) -> tuple[str, float]:
    """
    “不闭合时把差额最大的那项挑出来”的正式挑选规则：

    先拿三个有独立直接计量的科目 N/T/C 做账实核对；若任一科目的差额超过
    精确复核容差，挑其中差额最大的科目。只有三个科目都账实相符时，才把
    反解残差 R 本身（积分/未建模效应）作为最大差额项。

    不能直接把 |R| 与 δN/δT/δC 放一起取 max：R 是吸收项，某科目的记账错误
    会机械地搬进 R，两边几乎相等，直接 max 会永远偏向“残差”。
    """
    d = term_discrepancies(run)
    tol = exact_frac * run["E0"]
    channels = {k: d[k] for k in ("N", "T", "C")}
    key = max(channels, key=channels.get)
    if channels[key] > tol:
        return key, channels[key]
    return "R", abs(run["R"])


def canonical_impact_demo() -> Dict:
    """
    solve.c 自检同款单次碰撞：v=(vx=4, vy=-1)，n=+y，m=1，e=0，mu=1。
    给出三种记账口径，精确展示 “R≡0 且 T≡0 => 切向冲量没参与能量结算”。
    """
    m, vn, vt, e, mu = 1.0, -1.0, 4.0, 0.0, 1.0
    vn_post = -e * vn
    red = min(mu * abs(vn), abs(vt))
    vt_post = max(abs(vt) - red, 0.0)
    Ln = 0.5 * m * (vn * vn - vn_post * vn_post)   # 0.5 J
    Lt = 0.5 * m * (vt * vt - vt_post * vt_post)   # 3.5 J
    dE_full = -Ln - Lt                              # -4.0 J
    dE_vert = 0.5 * m * vn_post ** 2 - 0.5 * m * vn ** 2  # -0.5 J，接触点 PE 不变

    def book(obs_dE, T_book):
        N_book = -Ln
        C_book = 0.0
        R = obs_dE - N_book - T_book - C_book
        return {"N": N_book, "T": T_book, "C": C_book, "R": R,
                "sum4": N_book + T_book + C_book + R, "dE": obs_dE}

    return {
        "input": {"m": m, "vn_pre": vn, "vt_pre": vt, "e": e, "mu": mu,
                  "vn_post": vn_post, "vt_post": vt_post},
        "Ln_J": Ln, "Lt_J": Lt,
        "A_full_complete": book(dE_full, -Lt),       # 完整口径 + 摩擦入账
        "B_full_omitT": book(dE_full, 0.0),          # 全能量但漏摩擦账
        "C_vertical_omitT": book(dE_vert, 0.0),      # 截断竖直能量 + 漏摩擦账
    }


# --------------------------------------------------------------------------- #
# 粒度实验
# --------------------------------------------------------------------------- #
# 五档粒度（由粗到细，视图能力嵌套）
TIERS = [
    ("G0", "全程总账"),
    ("G1", "下落/碰撞/反弹 阶段账"),
    ("G2", "逐帧四项 + 速度系数公式复核"),
    ("G3", "帧内分段计量 (积分/冲量/修正)"),
    ("G4", "分量能量 + 冲量/动量交叉校验"),
]

# 每档每帧计量操作数的显式成本模型（能量求值/公式/不变量的标量运算加权）
TIER_COST_OPS = {
    "G0": 0.05,    # 全程 2 次能量求值摊到约 40 帧
    "G1": 0.15,    # 4 个阶段边界
    "G2": 6.0,     # 每帧 2 次 E + 接触帧 N/T 公式与状态差复核
    "G3": 9.0,     # 再加 3 个分段能量快照
    "G4": 13.0,    # 再加分量能、Jn/Jt 动量与摩擦锥不变量
}

CLOSE_TOL_PCT = 1.0  # “非闭行”的重要性阈值：1% E0


def _phases(run: Dict) -> Dict[str, List[Dict]]:
    cf = run["contact_frame"]
    fr = run["frames"]
    if cf < 0:
        return {"fall": fr, "contact": [], "rise": []}
    contact = [f for f in fr if f["contact"]]
    fall = [f for f in fr if f["k"] < cf]
    last_contact_k = max((f["k"] for f in contact), default=cf)
    rise = [f for f in fr if f["k"] > last_contact_k]
    return {"fall": fall, "contact": contact, "rise": rise}


def _sum(rows: List[Dict], key: str) -> float:
    return float(sum(r[key] for r in rows))


def nonclosing_stats(run: Dict) -> Dict[str, Dict[str, int]]:
    """
    各粒度视图下超过 1% E0 重要性闸的未闭合账行。

    rows   : 未闭合账行（阶段/帧/分段）数；
    frames : 这些账行覆盖的“原子仿真帧”数，即“非闭帧数”。粗粒度一行不过闸，
             该行覆盖的全部帧都无法被单独证明闭合（爆炸半径）；细粒度收敛到
             真正肇事的那一帧。G4 的能量行与 G3 相同，新增的是不变量检查。
    """
    tol = CLOSE_TOL_PCT / 100.0 * run["E0"]
    fr = run["frames"]
    ph = _phases(run)
    out: Dict[str, Dict[str, int]] = {}

    g0_bad = abs(run["R"]) > tol
    out["G0"] = {"rows": int(g0_bad), "frames": run["n_frames"] if g0_bad else 0}

    g1_rows = g1_frames = 0
    for rows in ph.values():
        if rows and abs(_sum(rows, "R")) > tol:
            g1_rows += 1
            g1_frames += len(rows)
    out["G1"] = {"rows": g1_rows, "frames": g1_frames}

    bad2 = [f for f in fr if abs(f["R"]) > tol]
    out["G2"] = {"rows": len(bad2), "frames": len(bad2)}

    g3_rows = 0
    bad3_frames = set()
    for f in fr:
        stage_bad = False
        if abs(f["dE_int"]) > tol:
            g3_rows += 1
            stage_bad = True
        if f["contact"]:
            r_solve = f["dE_solve"] - f["N_book"] - f["T_book"]
            if abs(r_solve) > tol:
                g3_rows += 1
                stage_bad = True
        if f["pen"] > 0.0 or abs(f["dE_corr"]) > 0.0:
            r_corr = f["dE_corr"] - f["C_book"]
            if abs(r_corr) > tol:
                g3_rows += 1
                stage_bad = True
        if stage_bad:
            bad3_frames.add(f["k"])
    out["G3"] = {"rows": g3_rows, "frames": len(bad3_frames)}
    out["G4"] = {"rows": g3_rows, "frames": len(bad3_frames)}
    return out


def count_nonclosing(run: Dict) -> Dict[str, int]:
    """向后兼容：只取未闭合行数。"""
    stats = nonclosing_stats(run)
    return {g: stats[g]["rows"] for g, _ in TIERS}


def tier_flags(run: Dict) -> Dict[str, set]:
    """
    每档粒度“可见的检查”会点名哪些账项 {N,T,C,R}。
    完全由真实记录与该档保留的字段决定，与故障编号无关；细档能力是粗档超集。
    一个缺陷算“可定位”，当且仅当它的目标账项出现在该档 flags 里。
    """
    fr = run["frames"]
    ph = _phases(run)
    E0 = run["E0"]
    exact = 1e-7 * E0   # 恒等式/公式复核容差（精确账应到浮点量级）
    mat = 0.01 * E0     # 1% E0 材料性闸
    fl: Dict[str, set] = {g: set() for g, _ in TIERS}

    # G0：只有一个总残差。再大也只知道“总账异常”，无法点名 N/T/C/R 科目。

    # G1：下落/反弹飞行阶段的残差累积越 1% 闸 -> 可点名 R（积分截断/未建模）
    flight_R = _sum(ph["fall"], "R") + _sum(ph["rise"], "R")
    if abs(flight_R) > mat:
        fl["G1"].add("R")

    # G2：逐帧四项，并用“账本自己保留的速度分量”做系数公式独立复核
    tangent_logged = not bool(run.get("flags", {}).get("trunc_E", False))
    for f in fr:
        if f["contact"]:
            if abs(f["N_book"] - f["N_true"]) > exact:
                fl["G2"].add("N")
            # 截断竖直口径的账本根本没保留 vt，这一项在 G2 无法复核 —— 正是 D3
            if tangent_logged and abs(f["T_book"] - f["T_true"]) > exact:
                fl["G2"].add("T")
        if f["pen"] > 0.0:
            c_expected = Params.m * Params.g * f["dy_logged"] + f.get("dK_clamp", 0.0)
            if abs(f["C_book"] - c_expected) > exact:
                fl["G2"].add("C")
        if abs(f["R"]) > mat:
            fl["G2"].add("R")
    fl["G2"] |= fl["G1"]   # 逐帧视图总能重新聚合出阶段视图

    # G3：帧内分段能量快照，能抓住“系数表式没覆盖的段内操作”
    fl["G3"] |= fl["G2"]
    for f in fr:
        if f["pen"] > 0.0 and abs(f["dE_corr"] - f["C_book"]) > exact:
            fl["G3"].add("C")
        if abs(f["dE_int"]) > exact:
            fl["G3"].add("R")   # 半隐式欧拉的分段截断误差显形

    # G4：独立物理状态/冲量交叉校验，不受账本能量口径（如漏 vx²）影响
    fl["G4"] |= fl["G3"]
    mu = run["params"]["mu"]
    for f in fr:
        if not f["contact"]:
            continue
        tangential_action = abs(f["vt_pre"] - f["vx"]) > 1e-9
        if tangential_action and abs(f["T_book"]) < exact:
            fl["G4"].add("T")   # Δvt≠0 却没有摩擦能量账 -> 切向冲量未结算
        cone_violation = abs(f["Jt"]) > mu * abs(f["Jn"]) + 1e-9
        if f["reversed_tangent"] or cone_violation:
            fl["G4"].add("T")   # 摩擦不反向 / |Jt|≤μ|Jn| 不变量
    return fl


def locate_fault(target: str, run: Dict) -> Dict[str, bool]:
    """保留旧入口：目标账项是否被各档粒度点名。"""
    fl = tier_flags(run)
    return {g: target in fl[g] for g, _ in TIERS}


# 8 个注入缺陷的真实运行配置
def fault_suite() -> List[Dict]:
    base = dict(H=4.0, vx0=4.0)
    p60 = Params(dt=1 / 60, mu=0.3)
    p30 = Params(dt=1 / 30, mu=0.3)
    p60grip = Params(dt=1 / 60, mu=0.8)
    return [
        {"id": "D1", "target": "N", "name": "法向损失把速度差 |vn|-|vn'| 当能量",
         "run": simulate(p=p60, **base, speed_N=True)},
        {"id": "D2", "target": "T", "name": "切向冲量已施加，摩擦能账记 0（全能量口径）",
         "run": simulate(p=p60, **base, no_T=True)},
        {"id": "D3", "target": "T", "name": "能量口径漏算 vx²，且摩擦账记 0（R≡0,T≡0）",
         "run": simulate(p=p60, **base, no_T=True, trunc_E=True)},
        {"id": "D4", "target": "C", "name": "位置修正注入的势能完全不入账",
         "run": simulate(p=p60, **base, no_C=True)},
        {"id": "D5", "target": "R", "name": "积分步长过粗（1/30 s），飞行段能量漂移累积",
         "run": simulate(p=p30, **base)},
        {"id": "D6", "target": "N", "name": "账本恢复系数用错（账本 e=0.9，实际 e=0.5）",
         "run": simulate(p=p60, **base, e_book=0.9)},
        {"id": "D7", "target": "C", "name": "位置修正实际执行两次，账只记一次",
         "run": simulate(p=p60, **base, double_C=True)},
        {"id": "D8", "target": "T", "name": "摩擦冲量不封项导致切向反向（能量可能仍闭合）",
         "run": simulate(p=p60grip, **base, reverse_friction=True)},
    ]


# --------------------------------------------------------------------------- #
# 主流程：四档落高表 + 粒度实验 + 落盘
# --------------------------------------------------------------------------- #
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="artifacts/energy-ledger")
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    heights = [1.0, 2.0, 4.0, 8.0]
    rows = []
    runs = {}
    for H in heights:
        r = simulate(H=H, vx0=0.0, p=Params(e=0.5, mu=0.2, dt=1e-3))
        runs[H] = r
        cm = closure_metrics(r)
        wk, wv = worst_discrepancy(r)
        rows.append({
            "H_m": H,
            "frames": r["n_frames"],
            "E0_J": r["E0"],
            "E_apex_J": r["E_apex_full"],
            "N_J": r["N"], "T_J": r["T"], "C_J": r["C"], "R_J": r["R"],
            "dE_J": r["dE"],
            "sum4_J": r["N"] + r["T"] + r["C"] + r["R"],
            "identity_err_pctE0": cm["identity_err_pctE0"],
            "residual_pctE0": cm["residual_pctE0"],
            "pass_1pct": cm["pass_1pct"],
            "worst_term": TERM_CN[wk], "worst_discrepancy_J": wv,
            "apex_height_m": r["z_apex"],
        })

    csv_path = os.path.join(args.out, "four_height_ledger.csv")
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for row in rows:
            w.writerow(row)

    # JSON（不含庞大 frames 的精简版 + 完整一份供复核）
    slim = []
    for H, r in runs.items():
        cm = closure_metrics(r)
        slim.append({k: v for k, v in r.items() if k != "frames"} | {"closure": cm})
    with open(os.path.join(args.out, "four_height_ledger.json"), "w",
              encoding="utf-8") as fh:
        json.dump(slim, fh, ensure_ascii=False, indent=2)

    # ---- 粒度实验 ----
    faults = fault_suite()
    loc_matrix: List[Dict] = []
    nonclose_rows_total = {g: 0 for g, _ in TIERS}
    nonclose_frames_total = {g: 0 for g, _ in TIERS}

    # 非闭帧计数包含四档标称落高 + 8 个故障
    for H in heights:
        st = nonclosing_stats(runs[H])
        for g in nonclose_rows_total:
            nonclose_rows_total[g] += st[g]["rows"]
            nonclose_frames_total[g] += st[g]["frames"]

    picker_rows = []
    for fdef in faults:
        r = fdef["run"]
        loc = locate_fault(fdef["target"], r)
        st = nonclosing_stats(r)
        for g in nonclose_rows_total:
            nonclose_rows_total[g] += st[g]["rows"]
            nonclose_frames_total[g] += st[g]["frames"]
        wk, wv = worst_discrepancy(r)
        d = term_discrepancies(r)
        picker_rows.append({
            "id": fdef["id"], "target": fdef["target"],
            "name": fdef["name"],
            "worst_term": wk, "worst_term_cn": TERM_CN[wk],
            "worst_discrepancy_J": wv,
            "dN_J": d["N"], "dT_J": d["T"], "dC_J": d["C"], "dR_J": d["R"],
            "R_book_J": r["R"], "T_book_J": r["T"],
            "R_fullstate_J": r["R_fullstate"],
            "E0_J": r["E0"],
            "worst_pctE0": 100 * wv / r["E0"],
            **{f"loc_{g}": int(loc[g]) for g, _ in TIERS},
            **{f"ncrows_{g}": st[g]["rows"] for g, _ in TIERS},
            **{f"ncframes_{g}": st[g]["frames"] for g, _ in TIERS},
        })
        loc_matrix.append({"id": fdef["id"], "target": fdef["target"],
                           "name": fdef["name"], **loc})

    located_count = {g: sum(1 for fm in loc_matrix if fm[g]) for g, _ in TIERS}

    gran_rows = []
    for g, label in TIERS:
        gran_rows.append({
            "tier": g, "label": label,
            "located_problems": located_count[g],
            "nonclosing_rows": nonclose_rows_total[g],
            "nonclosing_frames": nonclose_frames_total[g],
            "ops_per_frame": TIER_COST_OPS[g],
        })
    with open(os.path.join(args.out, "granularity_curve.csv"), "w",
              newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(gran_rows[0].keys()))
        w.writeheader()
        w.writerows(gran_rows)

    with open(os.path.join(args.out, "fault_picker.csv"), "w",
              newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(picker_rows[0].keys()))
        w.writeheader()
        w.writerows(picker_rows)

    demo = canonical_impact_demo()
    with open(os.path.join(args.out, "canonical_impact_demo.json"), "w",
              encoding="utf-8") as fh:
        json.dump(demo, fh, ensure_ascii=False, indent=2)

    with open(os.path.join(args.out, "localization_matrix.json"), "w",
              encoding="utf-8") as fh:
        json.dump({"tiers": TIERS, "matrix": loc_matrix,
                   "located_count": located_count,
                   "nonclosing_rows_total": nonclose_rows_total,
                   "nonclosing_frames_total": nonclose_frames_total},
                  fh, ensure_ascii=False, indent=2)

    # ---- dt 收敛扫描（辅助决定参考步长） ----
    dt_rows = []
    for H in heights:
        for dt in (1 / 30, 1 / 60, 1 / 120, 1 / 240, 1 / 480, 1 / 1000):
            r = simulate(H=H, p=Params(e=0.5, mu=0.2, dt=dt))
            cm = closure_metrics(r)
            dt_rows.append({"H_m": H, "dt_s": dt,
                            "residual_pctE0": cm["residual_pctE0"],
                            "pass_1pct": cm["pass_1pct"]})
    with open(os.path.join(args.out, "dt_sweep.csv"), "w",
              newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(dt_rows[0].keys()))
        w.writeheader()
        w.writerows(dt_rows)

    # ---- 控制台摘要 ----
    print("== 四档落高能量账 (J, m=1 kg, e=0.5, mu=0.2, dt=1e-3) ==")
    hdr = f"{'H':>4} {'E0':>8} {'Eapex':>8} {'N':>8} {'T':>6} {'C':>7} {'R':>8} {'|R|/E0%':>8} {'id%E0':>9} pass"
    print(hdr)
    for row in rows:
        print(f"{row['H_m']:>4.0f} {row['E0_J']:>8.3f} {row['E_apex_J']:>8.3f} "
              f"{row['N_J']:>8.3f} {row['T_J']:>6.3f} {row['C_J']:>7.4f} "
              f"{row['R_J']:>8.4f} {row['residual_pctE0']:>8.4f} "
              f"{row['identity_err_pctE0']:>9.2e} {int(row['pass_1pct'])}")

    print("\n== 粒度曲线 ==")
    for gr in gran_rows:
        print(f"{gr['tier']} {gr['label']:<32} locate={gr['located_problems']}/8 "
              f"nc_rows={gr['nonclosing_rows']:>3} nc_frames={gr['nonclosing_frames']:>4} "
              f"ops/f={gr['ops_per_frame']}")

    d = canonical_impact_demo()
    print("\n== 单次碰撞三口径 (vn=-1, vt=4, e=0, mu=1; Ln=0.5J Lt=3.5J) ==")
    for name in ("A_full_complete", "B_full_omitT", "C_vertical_omitT"):
        b = d[name]
        print(f"{name:<20} N={b['N']:+.2f} T={b['T']:+.2f} C={b['C']:+.2f} "
              f"R={b['R']:+.2e} sum4-dE={b['sum4'] - b['dE']:+.1e}")

    print("\n== 故障最大差额项 ==")
    for pr in picker_rows:
        print(f"{pr['id']} target={TERM_CN[pr['target']]:<5} pick={pr['worst_term_cn']:<5} "
              f"maxd={pr['worst_discrepancy_J']:8.3f}J ({pr['worst_pctE0']:6.2f}% E0) "
              f"R_book={pr['R_book_J']:8.3f} R_full={pr['R_fullstate_J']:8.4f}")

    if not args.no_plots:
        make_plots(args.out, rows, gran_rows, picker_rows, runs, demo)

    # 四档必须全部通过 1%
    return 0 if all(row["pass_1pct"] for row in rows) else 2


# --------------------------------------------------------------------------- #
# 绘图
# --------------------------------------------------------------------------- #
def make_plots(out_dir: str, rows: List[Dict], gran_rows: List[Dict],
               picker_rows: List[Dict], runs: Dict, demo: Dict) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager

    for fam in ("Noto Sans CJK SC", "WenQuanYi Zen Hei", "sans-serif"):
        try:
            font_manager.findfont(fam, fallback_to_default=False)
            plt.rcParams["font.sans-serif"] = [fam]
            break
        except Exception:
            continue
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["font.size"] = 10

    # ---------------- 图 1：四档能量瀑布 ----------------
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.6))
    axes = axes.ravel()
    for ax, row in zip(axes, rows):
        H = row["H_m"]
        comps = [
            ("初始势能", row["E0_J"], "total0"),
            ("法向恢复", row["N_J"], "N"),
            ("切向摩擦", row["T_J"], "T"),
            ("数值修正", row["C_J"], "C"),
            ("残差", row["R_J"], "R"),
            ("弹起机械能", row["E_apex_J"], "total1"),
        ]
        xs = list(range(len(comps)))
        cum = 0.0
        levels = []  # 每根柱的 (bottom, top, label_y)
        prev_top = None
        for i, (label, val, kind) in enumerate(comps):
            if kind.startswith("total"):
                ax.bar(i, val, bottom=0.0, color=INK2 if kind == "total0" else INK,
                       width=0.62, zorder=3)
                levels.append((0.0, val, val))
                cum = val
            else:
                new_cum = cum + val
                lo, hi = min(cum, new_cum), max(cum, new_cum)
                ax.bar(i, hi - lo, bottom=lo, color=TERM_COLOR[kind],
                       width=0.62, zorder=3)
                # 与上一柱末端的连接虚线
                ax.plot([i - 1 + 0.31, i - 0.31], [cum, cum],
                        color=MUTED, lw=0.8, ls="--", zorder=2)
                top = hi
                levels.append((lo, hi, top))
                cum = new_cum
        ax.set_xticks(xs)
        ax.set_xticklabels([c[0] for c in comps], rotation=20, ha="right", fontsize=8.5)
        id_err = row["identity_err_pctE0"]
        id_txt = "≈0（反解恒等）" if id_err == 0.0 else f"{id_err:.1e}%$E_0$"
        ax.set_title(f"落高 H = {H:.0f} m　|R|/$E_0$ = {row['residual_pctE0']:.3f}%　"
                     f"四项恒等误差 {id_txt}",
                     fontsize=10.5)
        ax.grid(axis="y", color=GRID, lw=0.7, zorder=0)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.spines["left"].set_color("#c3c2b7")
        ax.spines["bottom"].set_color("#c3c2b7")
        ax.tick_params(colors=INK2, labelsize=8.5)
        ax.axhline(0, color="#c3c2b7", lw=1.0)
        # 数值直标：总量柱标在柱顶上方，增量柱标在变化端外侧
        yspan = max(c[1] for c in comps) + 1.0
        for i, (label, val, kind) in enumerate(comps):
            lo, hi, top = levels[i]
            if kind.startswith("total"):
                ax.annotate(f"{val:.2f}", (i, hi), textcoords="offset points",
                            xytext=(0, 4), ha="center", fontsize=8, color=INK)
            else:
                signed_end = lo if val < 0 else hi
                off = -11 if val < 0 else 4
                ax.annotate(f"{val:+.3f}", (i, signed_end),
                            textcoords="offset points", xytext=(0, off),
                            ha="center", fontsize=7.6, color=INK2)
        ax.set_ylim(top=yspan * 1.14)
    fig.suptitle("四档落高能量账（m=1 kg，e=0.5，μ=0.2，dt=1 ms；单位 J）",
                 fontsize=13, color=INK, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    p1 = os.path.join(out_dir, "fig_four_height_waterfall.png")
    fig.savefig(p1, dpi=170, facecolor=SURFACE)
    plt.close(fig)

    # ---------------- 图 2：粒度两条曲线 + 开销 ----------------
    tiers = [g["tier"] for g in gran_rows]
    x = list(range(len(tiers)))
    loc = [g["located_problems"] for g in gran_rows]
    nc = [g["nonclosing_frames"] for g in gran_rows]
    nr = [g["nonclosing_rows"] for g in gran_rows]
    ops = [g["ops_per_frame"] for g in gran_rows]
    labels = [g["label"] for g in gran_rows]

    fig, (ax1, ax2, ax3) = plt.subplots(
        3, 1, figsize=(9.2, 9.4), sharex=True,
        gridspec_kw={"height_ratios": [1.15, 1.0, 0.8]})

    ax1.plot(x, loc, "-o", color="#2a78d6", lw=2, markersize=8, zorder=3)
    for xi, yi in zip(x, loc):
        ax1.annotate(f"{yi}/8", (xi, yi), textcoords="offset points",
                     xytext=(0, 9), ha="center", fontsize=10, color=INK)
    ax1.set_ylabel("可定位的注入问题数")
    ax1.set_ylim(-0.4, 8.8)
    ax1.set_yticks(range(0, 9, 2))
    ax1.set_title("记账粒度：定位能力 / 非闭行数 / 每帧开销（五档，由粗到细）",
                  fontsize=12.5)
    ax1.grid(axis="y", color=GRID, lw=0.7)
    ax1.axvline(2.0, color=MUTED, ls=":", lw=1.2)
    ax1.annotate("日常定 G2 + 常开切向冲量不变量 → 7/8\n"
                 "定位修正类(D7)时再开 G3 分段快照",
                 xy=(2, 5), xytext=(2.32, 0.7),
                 fontsize=9.3, color=INK2,
                 arrowprops=dict(arrowstyle="->", color=MUTED, lw=1.1))

    ax2.plot(x, nc, "-s", color="#eb6834", lw=2, markersize=7.5, zorder=3)
    for xi, yi, ri in zip(x, nc, nr):
        ax2.annotate(f"{yi} 帧\n({ri} 行)", (xi, yi), textcoords="offset points",
                     xytext=(0, 9 if yi < max(nc) * 0.8 else -26),
                     ha="center", fontsize=8.8, color=INK)
    ax2.set_ylabel("非闭帧数（未闭合账行覆盖帧）")
    ax2.set_ylim(bottom=-max(nc) * 0.06, top=max(nc) * 1.18 if max(nc) else 1)
    ax2.grid(axis="y", color=GRID, lw=0.7)

    ax3.plot(x, ops, "-^", color="#1baf7a", lw=2, markersize=8, zorder=3)
    for xi, yi in zip(x, ops):
        ax3.annotate(f"{yi:g}", (xi, yi), textcoords="offset points",
                     xytext=(0, 8), ha="center", fontsize=9.5, color=INK)
    ax3.set_ylabel("计量操作 / 帧")
    ax3.set_xticks(x)
    ax3.set_xticklabels([f"{t}\n{lab}" for t, lab in zip(tiers, labels)],
                        fontsize=8.6)
    ax3.grid(axis="y", color=GRID, lw=0.7)

    for ax in (ax1, ax2, ax3):
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.spines["left"].set_color("#c3c2b7")
        ax.spines["bottom"].set_color("#c3c2b7")
        ax.tick_params(colors=INK2)
    fig.tight_layout()
    p2 = os.path.join(out_dir, "fig_granularity_curves.png")
    fig.savefig(p2, dpi=170, facecolor=SURFACE)
    plt.close(fig)

    # ---------------- 图 3：切向冲量“隐身”诊断（标准单次碰撞三口径） ----------------
    groups = [
        ("口径 A\n全能量 + 摩擦入账", demo["A_full_complete"]),
        ("口径 B\n全能量但漏摩擦账", demo["B_full_omitT"]),
        ("口径 C\n只算竖直能量且漏摩擦账", demo["C_vertical_omitT"]),
    ]
    terms = ["N", "T", "C", "R"]
    fig, ax = plt.subplots(figsize=(9.4, 5.4))
    ng, nt = len(groups), len(terms)
    bw = 0.18
    for gi, (gname, book) in enumerate(groups):
        center = gi * (nt * bw + 0.9)
        for ti, tk in enumerate(terms):
            xpos = center + (ti - 1.5) * bw
            val = book[tk]
            ax.bar(xpos, val, width=bw * 0.92, color=TERM_COLOR[tk], zorder=3,
                   label=TERM_CN[tk] if gi == 0 else None)
            if abs(val) > 1e-9:
                off = 4 if val >= 0 else -12
                ax.annotate(f"{val:+.2f}", (xpos, val), textcoords="offset points",
                            xytext=(0, off), ha="center", fontsize=8.6, color=INK2)
            else:
                ax.annotate("0", (xpos, 0), textcoords="offset points",
                            xytext=(0, 3), ha="center", fontsize=8.0, color=MUTED)
        dE = book["dE"]
        ax.annotate(f"总变化 ΔE={dE:+.2f}", (center, -4.35), ha="center",
                    fontsize=9.2, color=INK,
                    fontweight="bold" if gi == 2 else "normal")
        ax.text(center, 0.62, gname, ha="center", va="bottom", fontsize=9.3,
                color=INK)
    ax.axhline(0, color="#c3c2b7", lw=1.1)
    ax.set_xticks([])
    ax.set_ylabel("各项能量变化贡献 (J，损失为负)")
    ax.set_ylim(-4.9, 1.15)
    ax.set_title("同一个碰撞 vx: 4→3 m/s、vy: −1→0（e=0, μ=1）：三种记账口径\n"
                 "口径 C 里 R≡0 且 T≡0 —— 不是没有摩擦，是切向自由度根本没进能量账",
                 fontsize=12)
    ax.grid(axis="y", color=GRID, lw=0.7, zorder=0)
    for sp in ("top", "right", "bottom"):
        ax.spines[sp].set_visible(False)
    ax.spines["left"].set_color("#c3c2b7")
    ax.tick_params(colors=INK2)
    ax.legend(frameon=False, fontsize=9.5, ncol=4, loc="upper center",
              bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout()
    p3 = os.path.join(out_dir, "fig_tangent_diagnostic.png")
    fig.savefig(p3, dpi=170, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)

    print(f"\nplots: {p1}\n       {p2}\n       {p3}")


if __name__ == "__main__":
    raise SystemExit(main())
