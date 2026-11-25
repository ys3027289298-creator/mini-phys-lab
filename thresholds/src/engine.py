"""
engine.py — 把 collide / solve / integrate 三个 C 模块移植成 2D 仿真，
外加休眠/退化/穿透三套判定、带噪声观测与真值标签。

物理与 C 源码一一对应:
  collide.c : sphere_vs_aabb()  最近点、法线、深度、质心入盒走最浅轴
  solve.c   : resolve_contact() 法向冲量 + 切向摩擦，摩擦不反向
  integrate : step_body()       半隐式欧拉 + 位置修正，撞地清零法向速
所有观测值 = 真值 + 高斯噪声；误判统计始终用无噪真值打标签。
"""
import math
import random


# ---------------------------------------------------------------- 基础类型
class V2:
    __slots__ = ("x", "y")

    def __init__(self, x=0.0, y=0.0):
        self.x, self.y = x, y

    def __repr__(self):
        return f"V2({self.x:.4g},{self.y:.4g})"


class Aabb:
    def __init__(self, x0, y0, x1, y1):
        self.min = V2(x0, y0)
        self.max = V2(x1, y1)


class Ball:
    def __init__(self, x, y, r, restitution=0.3, friction=0.4):
        self.x, self.y, self.r = x, y, r
        self.vx = self.vy = 0.0
        self.px, self.py = x, y
        self.e = restitution
        self.mu = friction
        self.asleep = False
        self.sleep_timer = 0          # 连续低速帧计数
        self.wake_latency = 0         # >0 表示已被撞、等待唤醒的剩余帧
        self.mass = r * r             # 密度均匀，质量 ~ r^2

    @property
    def speed(self):
        return math.hypot(self.vx, self.vy)


# ---------------------------------------------------------------- collide.c 移植
def sphere_vs_aabb(b: Ball, box: Aabb):
    """返回 (hit, nx, ny, depth)；质心入盒时沿最浅轴推出，同 collide.c。"""
    qx = min(max(b.x, box.min.x), box.max.x)
    qy = min(max(b.y, box.min.y), box.max.y)
    dx, dy = b.x - qx, b.y - qy
    d2 = dx * dx + dy * dy
    if d2 > b.r * b.r:
        return 0, 0.0, 0.0, 0.0
    if d2 > 1e-12:
        d = math.sqrt(d2)
        return 1, dx / d, dy / d, b.r - d
    px = min(b.x - box.min.x, box.max.x - b.x)
    py = min(b.y - box.min.y, box.max.y - b.y)
    if px <= py:
        nx = -1.0 if b.x < 0.5 * (box.min.x + box.max.x) else 1.0
        return 1, nx, 0.0, px + b.r
    ny = -1.0 if b.y < 0.5 * (box.min.y + box.max.y) else 1.0
    return 1, 0.0, ny, py + b.r


def sphere_vs_sphere(a: Ball, b: Ball):
    dx, dy = b.x - a.x, b.y - a.y
    d = math.hypot(dx, dy)
    rsum = a.r + b.r
    if d >= rsum or d < 1e-12:
        if d < 1e-12:
            return 1, 1.0, 0.0, rsum          # 完全重叠，退化到 x 轴
        return 0, 0.0, 0.0, 0.0
    return 1, dx / d, dy / d, rsum - d


# ---------------------------------------------------------------- solve.c 移植
def resolve_contact(va, vb, nx, ny, e, mu, inv_mass_sum):
    """
    对两物体（vb 可为 None 表示静世界）施加冲量，返回法向相对速率 |vn|。
    法向: new_vn = -e*vn；切向: 衰减 mu*|vn|，clamp 到不反向——与 solve.c 同式。
    """
    if vb is None:
        vbx = vby = 0.0
        ima = inv_mass_sum
        imb = 0.0
    else:
        vbx, vby = b_vx, b_vy = vb.vx, vb.vy
        ima, imb = 1.0 / va.mass, 1.0 / vb.mass
    rvx, rvy = va.vx - vbx, va.vy - vby
    vn = rvx * nx + rvy * ny
    if vn > 0.0:
        return 0.0                                  # 已分离
    tx, ty = rvx - vn * nx, rvy - vn * ny
    tl = math.hypot(tx, ty)
    new_vn = -vn * e
    reduce = min(mu * abs(vn), tl)
    scale = (tl - reduce) / tl if tl > 1e-12 else 0.0
    jx = tx * scale + nx * new_vn - rvx            # 相对速度的总改变量
    jy = ty * scale + ny * new_vn - rvy
    va.vx += jx * ima / (ima + imb)
    va.vy += jy * ima / (ima + imb)
    if vb is not None:
        vb.vx -= jx * imb / (ima + imb)
        vb.vy -= jy * imb / (ima + imb)
    return abs(vn)


# ---------------------------------------------------------------- 世界与判定
class Thresholds:
    def __init__(self, v_sleep, d_deg, p_tun, k_frames=4, wake_latency=0):
        self.v_sleep = v_sleep
        self.d_deg = d_deg
        self.p_tun = p_tun
        self.k_frames = k_frames
        self.wake_latency = wake_latency


class Counters:
    def __init__(self):
        # 休眠
        self.false_sleep = 0          # 真在动却睡着（帧）
        self.false_awake = 0          # 真静息却睡不着（静息片段）
        self.rest_episodes = 0
        self.moving_frames = 0        # 分母: true_v > v_min 的球-帧
        self.rest_frames = 0          # 分母: true_v < 0.5 v_min 的球-帧
        # 退化
        self.deg_miss = 0             # 真重叠 >= delta_res 却被退化
        self.deg_false = 0            # 真重叠 < delta_res 却上完整求解
        self.deep_probes = 0          # 分母: 真重叠 >= delta_res
        self.shallow_probes = 0       # 分母: 真重叠 < delta_res
        # 穿透
        self.pen_miss = 0             # 真风险未触发 CCD
        self.pen_false = 0            # 无风险却触发 CCD
        self.risk_probes = 0          # 分母: 真穿透风险
        self.safe_probes = 0          # 分母: 无风险接触
        self.ccd_triggers = 0
        # 耦合
        self.push_through = 0         # 休眠球被推穿（C1）
        self.jitter_alarms = 0        # 醒着的静息球抖动诱发误报（C2）
        self.never_slept = 0          # 慢速球整段没睡着（C3 片段）

    def as_dict(self):
        return {k: v for k, v in self.__dict__.items()}

    def merge(self, other):
        for k, v in other.as_dict().items():
            setattr(self, k, getattr(self, k) + v)

    def rates(self):
        """六个误判率，分母为 0 时记 None。误睡/误醒均为球-帧级，
        与 derive.py 中单帧概率经 k 帧组合的解析模型同口径。"""
        def safe(n, d):
            return None if d == 0 else n / d
        return {
            "false_sleep": safe(self.false_sleep, self.moving_frames),
            "false_awake": safe(self.false_awake,
                                max(self.rest_frames, 1)),
            "deg_miss": safe(self.deg_miss, self.deep_probes),
            "deg_false": safe(self.deg_false, self.shallow_probes),
            "pen_miss": safe(self.pen_miss, self.risk_probes),
            "pen_false": safe(self.pen_false, self.safe_probes),
        }


class World:
    def __init__(self, obs, thr: Thresholds, seed=7, slab_h=None):
        self.o = obs
        self.thr = thr
        self.rng = random.Random(seed)
        self.balls = []
        g = obs.g
        # 默认水平地板；坡面场景把重力旋转即可（接触代码仍是球-盒）
        self.gx, self.gy = 0.0, -g
        h = slab_h if slab_h is not None else obs.h
        self.slab = Aabb(-50.0, -h, 50.0, 0.0)      # 薄板，顶面 y=0
        self.cnt = Counters()
        self._rest_open = {}                         # ball -> 连续静息帧

    def gauss(self, sigma):
        return self.rng.gauss(0.0, sigma)

    def add_ball(self, x, y, r=None, **kw):
        b = Ball(x, y, r if r is not None else self.o.r, **kw)
        self.balls.append(b)
        self._rest_open[id(b)] = 0
        return b

    def set_slope(self, theta):
        """等效坡面: 重力旋转 theta，地板仍水平。"""
        g = self.o.g
        self.gx, self.gy = g * math.sin(theta), -g * math.cos(theta)

    # ---- 带噪观测 -------------------------------------------------------
    def obs_speed(self, b):
        return abs(b.speed + self.gauss(self.o.sigma_v))

    def obs_depth(self, depth, sigma):
        return max(0.0, depth + self.gauss(sigma))

    # ---- 单帧 -----------------------------------------------------------
    def step(self):
        o, thr, cnt = self.o, self.thr, self.cnt
        slab = self.slab

        for b in self.balls:
            # 休眠判定（观测速 + 连续 k 帧）
            v_obs = self.obs_speed(b)
            true_v = b.speed
            v_min = o.g * o.dt * o.mu_s / math.sqrt(1.0 + o.mu_s ** 2)
            in_contact, _, _, depth0 = sphere_vs_aabb(b, slab)
            resting = in_contact and depth0 <= o.delta0 * 3.0 and true_v < 0.5 * v_min

            if not b.asleep and b.wake_latency == 0:
                if v_obs < thr.v_sleep:
                    b.sleep_timer += 1
                else:
                    b.sleep_timer = 0
                if b.sleep_timer >= thr.k_frames:
                    b.asleep = True
                    b.vx = b.vy = 0.0

            # 真值标签：动 / 静息（静息须在接触态，自由落体 v≈0 不算）
            if true_v > v_min:
                cnt.moving_frames += 1
                if b.asleep:
                    cnt.false_sleep += 1
            if resting:
                cnt.rest_frames += 1
                self._rest_open[id(b)] += 1
                # 帧级误醒：物理上已静息，但判定还没让它睡着
                if not b.asleep and b.wake_latency == 0 \
                        and self._rest_open[id(b)] >= thr.k_frames:
                    cnt.false_awake += 1
            else:
                if self._rest_open[id(b)] >= thr.k_frames and not b.asleep \
                        and b.wake_latency == 0:
                    cnt.rest_episodes += 1
                self._rest_open[id(b)] = 0

            # integrate.c 移植：休眠体不积分（唤醒延迟内亦然）
            if not b.asleep and b.wake_latency == 0:
                b.px, b.py = b.x, b.y
                b.vx += self.gx * o.dt
                b.vy += self.gy * o.dt
                b.x += b.vx * o.dt
                b.y += b.vy * o.dt
            else:
                b.px, b.py = b.x, b.y

        # 球-球接触（collide + solve 移植），休眠球被撞 -> 登记唤醒
        for i in range(len(self.balls)):
            a = self.balls[i]
            for j in range(i + 1, len(self.balls)):
                b = self.balls[j]
                hit, nx, ny, depth = sphere_vs_sphere(a, b)
                if not hit:
                    continue
                if a.asleep or b.asleep:
                    sleeper = a if a.asleep else b
                    striker = b if a.asleep else a
                    if striker.speed > thr.v_sleep:
                        if sleeper.wake_latency == 0 and sleeper.asleep:
                            sleeper.wake_latency = thr.wake_latency + 1
                        sleeper.asleep = False
                # 退化判定：观测深度 < d_deg 只走法向冲量
                d_obs = self.obs_depth(depth, o.sigma_d)
                if depth >= o.eps_r * o.r:
                    cnt.deep_probes += 1
                else:
                    cnt.shallow_probes += 1
                if d_obs < thr.d_deg:
                    if depth >= o.eps_r * o.r:
                        cnt.deg_miss += 1
                    resolve_contact(a, b, nx, ny, a.e, 0.0,
                                    1.0 / a.mass + 1.0 / b.mass)
                else:
                    if depth < o.eps_r * o.r:
                        cnt.deg_false += 1
                    resolve_contact(a, b, nx, ny, a.e, min(a.mu, b.mu),
                                    1.0 / a.mass + 1.0 / b.mass)
                # 位置修正：按质量对半推开
                corr = depth
                wa = (1.0 / a.mass) / (1.0 / a.mass + 1.0 / b.mass)
                a.x -= nx * corr * wa
                a.y -= ny * corr * wa
                b.x += nx * corr * (1.0 - wa)
                b.y += ny * corr * (1.0 - wa)

        # 球-薄板接触（integrate.c 的地板修正 + 穿透 CCD）
        p_head = o.h - o.v_max * o.dt
        for b in self.balls:
            hit, nx, ny, depth = sphere_vs_aabb(b, slab)
            if not hit:
                # 整块越过薄板：本帧漏检的穿透（只可能在位移余量被吃干净时）
                if b.py - b.r > slab.max.y and b.y + b.r < slab.min.y:
                    cnt.pen_miss += 1
                    cnt.risk_probes += 1
                continue

            d_obs = self.obs_depth(depth, o.sigma_p)
            true_risk = depth > p_head
            in_slop = depth <= o.delta0            # 静息 slop，不属穿透判定域
            if true_risk:
                cnt.risk_probes += 1
            elif not in_slop:
                cnt.safe_probes += 1
            if d_obs > thr.p_tun:
                cnt.ccd_triggers += 1
                if not true_risk and not in_slop:
                    cnt.pen_false += 1
                if true_risk:
                    # CCD 扫掠钳制：回到本帧入射点
                    b.x, b.y = b.px, b.py
                    inward = b.vx * nx + b.vy * ny
                    if inward < 0.0:
                        b.vx -= inward * nx
                        b.vy -= inward * ny
            else:
                if true_risk:
                    cnt.pen_miss += 1

            # 耦合 C1：休眠/待唤醒球被位置修正推进薄板深处
            if (b.asleep or b.wake_latency > 0) and depth > p_head:
                cnt.push_through += 1
            # 耦合 C2：本该睡着却醒着的球，slop 上方抖动诱发 CCD 误报
            if not b.asleep and not in_slop and d_obs > thr.p_tun \
                    and not true_risk and b.speed < 0.5 * thr.v_sleep:
                cnt.jitter_alarms += 1

            # 静摩擦平衡：库仑上限 mu_s * 法向重力冲量，临界坡面上恰好等于切向冲量
            if not b.asleep and b.wake_latency == 0:
                tv = b.vx * nx + b.vy * ny
                tx, ty = b.vx - tv * nx, b.vy - tv * ny
                tl = math.hypot(tx, ty)
                normal_impulse = abs((self.gx * nx + self.gy * ny) * o.dt)
                hold = min(o.mu_s * normal_impulse, tl)
                if tl > 1e-12:
                    s = max(0.0, (tl - hold) / tl)
                    b.vx = tx * s + tv * nx
                    b.vy = ty * s + tv * ny

            # integrate.c 位置修正 + 法向速清零（slop 内保留，静息体因此稳定在 delta0）
            if depth > o.delta0:
                b.x += nx * (depth - o.delta0)
                b.y += ny * (depth - o.delta0)
            vn = b.vx * nx + b.vy * ny
            if vn < 0.0:
                b.vx -= vn * nx
                b.vy -= vn * ny

            if b.wake_latency > 0:
                b.wake_latency -= 1

    def run(self, frames):
        for _ in range(frames):
            self.step()
        # 收尾：仍开着的静息片段若最终没睡着，计入"永不休眠"耦合事件
        for b in self.balls:
            if self._rest_open[id(b)] >= self.thr.k_frames and not b.asleep:
                self.cnt.rest_episodes += 1
                self.cnt.never_slept += 1
        return self.cnt


# --------------------------------------------------------------- 受控判定探针
def contact_probes(obs, thr, depths, sigma, seed=11):
    """S3/S4 用：在指定真深度上注入带噪接触，直接统计阈值判定的混淆矩阵。"""
    rng = random.Random(seed)
    cnt = Counters()
    d_res = obs.eps_r * obs.r
    p_head = obs.h - obs.v_max * obs.dt
    for kind, true_d in depths:
        d_obs = max(0.0, true_d + rng.gauss(0.0, sigma))
        if kind == "degrade":
            if true_d >= d_res:
                cnt.deep_probes += 1
            else:
                cnt.shallow_probes += 1
            if d_obs < thr.d_deg:
                if true_d >= d_res:
                    cnt.deg_miss += 1
            elif true_d < d_res:
                cnt.deg_false += 1
        else:
            true_risk = true_d > p_head
            if true_risk:
                cnt.risk_probes += 1
            else:
                cnt.safe_probes += 1
            if d_obs > thr.p_tun:
                cnt.ccd_triggers += 1
                if not true_risk:
                    cnt.pen_false += 1
            elif true_risk:
                cnt.pen_miss += 1
    return cnt
