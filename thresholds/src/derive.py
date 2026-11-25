"""
derive.py — 三个阈值的闭式推导与目标误判率反解。

记号与可观测物理量（SI 单位）:
  r        球半径                     [m]      几何可观测量
  mu_s     最大静摩擦系数             [-]      材料可观测量
  g        重力加速度                 [m/s^2]
  dt       帧步长                     [s]
  sigma_v  观测速的 1σ 噪声           [m/s]    标定可观测量
  eps_r    最小可分辨重叠的相对量级   [-]      float32 余量 + 几何公差标定
  sigma_d  观测重叠深度的 1σ 噪声     [m]
  h        地板(薄板)厚度             [m]
  v_max    场景中单球最大速率         [m/s]
  sigma_p  观测穿透深度的 1σ 噪声     [m]
  delta0   静息接触的位置修正 slop    [m]

三个阈值共用同一个统计模型:
  物理尺度 s 由可观测量推出; 观测量 X ~ N(s, sigma^2);
  判定量 T 大于/小于 X 触发某类误判, 误判率随 T 单调,
  用正态分位数 z = Phi^-1(alpha) 反解阈值边界。
"""
from dataclasses import dataclass
import math


def norm_cdf(x: float) -> float:
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def norm_ppf(p: float) -> float:
    """Acklam 逆正态 CDF，全定义域相对误差 < 1.2e-7，足够阈值反解。"""
    if not 0.0 < p < 1.0:
        raise ValueError("norm_ppf needs 0 < p < 1")
    a = [-3.969683028665376e+01, 2.209460984245205e+02,
         -2.759285104469687e+02, 1.383577518672690e+02,
         -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02,
         -1.556989798598866e+02, 6.680131188771972e+01,
         -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01,
         -2.400758277161838e+00, -2.549732539343734e+00,
         4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01,
         2.445134137142996e+00, 3.754408661907416e+00]
    plow, phigh = 0.02425, 1.0 - 0.02425
    if p < plow:
        q = math.sqrt(-2.0 * math.log(p))
        x = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0)
    elif p <= phigh:
        q = p - 0.5
        r = q * q
        x = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / \
            (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1.0)
    else:
        q = math.sqrt(-2.0 * math.log(1.0 - p))
        x = -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
             ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0)
    # 一步 Halley 精修
    e = norm_cdf(x) - p
    u = e * math.sqrt(2.0 * math.pi) * math.exp(x*x/2.0)
    x -= u / (1.0 + x*u/2.0)
    return x


class _Norm:
    cdf = staticmethod(norm_cdf)
    ppf = staticmethod(norm_ppf)


norm = _Norm()


@dataclass
class Observables:
    r: float = 0.5            # 球半径 m（与 integrate 自测一致）
    mu_s: float = 0.4         # 静摩擦系数
    g: float = 9.8
    dt: float = 1.0 / 60.0
    sigma_v: float = 0.01     # 1 cm/s，滤波后的速率噪声
    eps_r: float = 2e-4       # 最小可分辨重叠 ≈ 半径的 2e-4（float32 余量之上的工程公差）
    sigma_d: float = 2e-5     # 0.02 mm 深度噪声
    h: float = 0.1            # 地板厚 10 cm
    v_max: float = 4.0        # 最大球速
    sigma_p: float = 3e-3     # 3 mm 穿透深度噪声
    delta0: float = 1e-3      # 静息 slop 1 mm
    k_frames: int = 4         # 连续 k 帧低速才休眠


# ---------------------------------------------------------------------------
# 1) 休眠阈值 v_sleep —— 从最大静摩擦能维持的坡度推
# ---------------------------------------------------------------------------
def limiting_slope(mu_s: float) -> float:
    """最大静摩擦能维持的坡度角 theta_max: tan(theta_max) = mu_s。"""
    return math.atan(mu_s)


def sleep_physical_scale(o: Observables) -> float:
    """
    临界坡面上，静摩擦恰好兜不住时球一帧内获得的最小物理速率:
        a_s   = g sin(theta_max) = g mu_s / sqrt(1 + mu_s^2)
        v_min = a_s * dt
    真在动的球，速率不可能长期低于 v_min；低于它的"运动"只可能是观测噪声。
    """
    sin_t = o.mu_s / math.sqrt(1.0 + o.mu_s ** 2)
    return o.g * o.dt * sin_t


def sleep_window(o: Observables, alpha_sleep: float, alpha_awake: float):
    """
    反解休眠阈值可行窗口。
    误睡（临界运动球被休眠）单帧率 p_fs = Phi((T - v_min)/sigma_v)，随 T 单调增 -> 上限:
        T <= v_min + sigma_v * Phi^-1(alpha_sleep)
    误醒（静息球睡不着），连续 k 帧判据:
        p_fa = 1 - Phi(T/sigma_v)^k  -> 下限 T >= sigma_v * Phi^-1((1-alpha_awake)^(1/k))
    返回 (下限, 上限, 物理尺度)；下限 > 上限 说明传感器/误判预算不可行。
    """
    v_min = sleep_physical_scale(o)
    upper = v_min + o.sigma_v * norm.ppf(alpha_sleep)
    lower = o.sigma_v * norm.ppf((1.0 - alpha_awake) ** (1.0 / o.k_frames))
    return lower, upper, v_min


def sleep_pick(o: Observables, alpha: float):
    """与实验校准同口径的反解取值：双侧等 alpha，取窗口 75% 分位。"""
    lo, hi, s = sleep_window(o, alpha, alpha)
    return lo + 0.75 * (hi - lo), lo, hi, s


def sleep_rates(o: Observables, T: float):
    """给定阈值的理论误睡/误醒率，供仿真对照。"""
    v_min = sleep_physical_scale(o)
    p_false_sleep = norm.cdf((T - v_min) / o.sigma_v)
    p_per_frame_awake = 1.0 - norm.cdf(T / o.sigma_v)
    p_false_awake = 1.0 - (1.0 - p_per_frame_awake) ** o.k_frames
    return p_false_sleep, p_false_awake


# ---------------------------------------------------------------------------
# 2) 退化阈值 d_deg —— 从半径与最小可分辨重叠推
# ---------------------------------------------------------------------------
def degrade_physical_scale(o: Observables) -> float:
    """
    必须保留完整流形求解的最小重叠:
        delta_res = eps_r * r
    低于它的接触在 float32 深度运算里不可分辨，退化为法向单冲量即可。
    """
    return o.eps_r * o.r


def degrade_window(o: Observables, alpha_miss: float, alpha_false: float):
    """
    漏判（真重叠 >= delta_res 却被退化）p = Phi((T - delta_res)/sigma_d)，随 T 单调增 -> 上限:
        T <= delta_res + sigma_d * Phi^-1(alpha_miss)
    误留（已分离/零重叠仍上完整求解），零重叠观测量 ~ N(0, sigma_d^2):
        下限 T >= sigma_d * Phi^-1(1 - alpha_false)
    """
    d_res = degrade_physical_scale(o)
    upper = d_res + o.sigma_d * norm.ppf(alpha_miss)
    lower = o.sigma_d * norm.ppf(1.0 - alpha_false)
    return lower, upper, d_res


def degrade_rates(o: Observables, T: float):
    d_res = degrade_physical_scale(o)
    p_miss = norm.cdf((T - d_res) / o.sigma_d)
    p_false = 1.0 - norm.cdf(T / o.sigma_d)
    return p_miss, p_false


# ---------------------------------------------------------------------------
# 3) 穿透阈值 p_tun —— 从地板厚度与单帧最大位移推
# ---------------------------------------------------------------------------
def penetration_physical_scale(o: Observables) -> float:
    """
    一帧最多走 delta_max = v_max * dt；要保证本帧不整块穿过厚 h 的地板，
    触发 CCD 的深度预算必须留下完整位移余量:
        p_head = h - v_max * dt
    """
    return o.h - o.v_max * o.dt


def penetration_window(o: Observables, alpha_miss: float, alpha_false: float,
                       coupling_margin: float = 0.0):
    """
    漏判（真实穿透风险未触发 CCD）随 T 单调增 -> 噪声反解上限:
        T <= p_head + sigma_p * Phi^-1(alpha_miss)
    误报（静息 slop 接触误触发昂贵 CCD），静息深度 ~ N(delta0, sigma_p^2):
        下限 T >= delta0 + sigma_p * Phi^-1(1 - alpha_false)
    耦合项: 休眠球被推、唤醒有延迟时，必须在独立预算之外额外留出
    coupling_margin（= 唤醒延迟帧数 * v_max * dt）的位移头裕，
    几何上限同步收窄；窗口倒置时 feasible=False 交上层定夺。
    """
    p_head = penetration_physical_scale(o)
    head = max(0.0, p_head - coupling_margin)
    upper_noise = head + o.sigma_p * norm.ppf(alpha_miss)
    upper = min(upper_noise, head)
    lower = o.delta0 + o.sigma_p * norm.ppf(1.0 - alpha_false)
    return lower, upper, p_head


def penetration_rates(o: Observables, T: float):
    p_head = penetration_physical_scale(o)
    p_miss = norm.cdf((T - p_head) / o.sigma_p)
    p_false = 1.0 - norm.cdf((T - o.delta0) / o.sigma_p)
    return p_miss, p_false


# ---------------------------------------------------------------------------
# 统一反解: 给定可接受误判率，定出阈值（上限优先，同时校验下限可行性）
# ---------------------------------------------------------------------------
def derive_thresholds(o: Observables, alpha: float,
                      wake_latency_frames: int = 0):
    """
    三类误判预算统一取 alpha（两侧等权）。返回每个阈值的取值、可行窗口与物理尺度。
    取值规则: 取窗口 75% 分位（在满足漏判上限的前提下尽量压低误报），
    并把结果 clip 到窗口内；窗口不可行时 feasible=False 交调用方处理。
    唤醒延迟按 v_max*dt*L 折算成穿透窗口的耦合位移裕量。
    """
    margin = o.v_max * o.dt * wake_latency_frames
    out = {}

    lo, hi, s = sleep_window(o, alpha, alpha)
    pick = lo + 0.75 * (hi - lo)
    out["sleep"] = dict(value=pick, lower=lo, upper=hi, scale=s, feasible=lo <= hi)

    lo, hi, s = degrade_window(o, alpha, alpha)
    pick = lo + 0.75 * (hi - lo)
    out["degrade"] = dict(value=pick, lower=lo, upper=hi, scale=s, feasible=lo <= hi)

    lo, hi, s = penetration_window(o, alpha, alpha, coupling_margin=margin)
    pick = lo + 0.75 * (hi - lo)
    out["penetrate"] = dict(value=max(0.0, pick), lower=lo, upper=hi,
                            scale=s, feasible=lo <= hi,
                            coupling_margin=margin)

    return out


if __name__ == "__main__":
    obs = Observables()
    print(f"theta_max          = {math.degrees(limiting_slope(obs.mu_s)):.2f} deg")
    print(f"v_min (sleep)      = {sleep_physical_scale(obs):.4f} m/s")
    print(f"delta_res (degrade)= {degrade_physical_scale(obs):.2e} m")
    print(f"p_head (penetrate) = {penetration_physical_scale(obs):.4f} m")
    for a in (0.05, 0.01, 0.005):
        print(f"\n--- alpha = {a:.3f} ---")
        for name, d in derive_thresholds(obs, a).items():
            flag = "OK " if d["feasible"] else "INF"
            print(f"  {name:9s} [{flag}] T={d['value']:.3e}  "
                  f"window=[{d['lower']:.3e}, {d['upper']:.3e}]  scale={d['scale']:.3e}")
