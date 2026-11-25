# collide_ss

球-球接触生成（法线 + 穿透深度），专门处理**退化法线**：球心距为零或小到
float32 无法分辨时，朴素 `n = Δp / |Δp|` 会 0/0。

## 构建与复现

```
make test    # 自检（开启 FE_DIVBYZERO/FE_INVALID，残留 0/0 直接 SIGFPE）
make repro   # 构建复现制品（不自动运行崩溃用例）
./naive_repro   # 朴素实现：预期 SIGFPE，exit=136
./nan_demo      # 默认 IEEE 语义：0/0 吐 NaN 并污染下游冲量
make plot    # 20 万样本阈值扫描 + 双曲线图（results/threshold_curves.png）
```

## 1. 崩溃复现与最小输入

朴素代码 `src/naive_repro.c`：

| 输入 | 结果 |
|---|---|
| 原始输入，原点附近 d=0.02 | float 可分辨，normal=(1,0,0)，**不崩** |
| 同一对球平移到 L=1e6 | 相减消去，d_flt=0，**0/0 → SIGFPE (exit=136)** |
| 最小输入：d=0、rA=rB=0.5 | 同因，0/0 |

原始输入在小坐标下没触发，但失效机理是"分母落到 0"——大坐标平移把同一对
球逼进该路径。不开浮点陷阱时（`nan_demo`），0/0 不报错而返回 NaN，下游
`solve` 的法向速度 `v·n`、单位法线自检全部被 NaN 静默污染（NaN 比较恒假，
自检防线也被绕过）。

## 2. 阈值判据（由半径与重叠深度反解）

等半径 r=0.5 时 rsum=1，穿透深度 δ = rsum − d = 1 − d。

- d < eps_t 时位置差无法确定法线方向；此时 δ > 1 − eps_t。
- 反解：**eps_t = rsum − δ_max**，即"仍允许用位置定法线"的最小球心距，
  等价于"能容忍多深的深穿透仍不走兜底"。
- 数值下界：float32 在坐标尺度 L 的分辨力 ULP(L) ≈ 1.2e-7·L。d 低于该
  尺度时减法消去，d 本身是 0 或噪声，eps_t 必须高于这条噪声地板。
  L=1 时地板 ≈1.2e-7；L=1e6 时地板抬到 ≈0.0625。

## 3. 三级兜底与反例（光调阈值不够）

`balls_contact` 的兜底链（`src/collide_ss.c`）：

1. **位置法线** `normalize(B.c − A.c)` —— 几何方向；
2. **相对速度法线** `normalize(A.v − B.v)` —— 位置退化时的运动学兜底；
3. **固定轴 (1,0,0)** —— 永不退化。

反例：两球球心完全重合**且速度相同**（自检 `counterexample_coincident_equal_velocity`）。
此时相对速度也是零向量，第 2 级同样 0/0——只调第 1 级阈值、或只做速度
兜底都救不了，必须落到固定轴。任何接触输出前过 `valid_unit` 不变量检查，
非有限值或非单位长度则拒绝输出（返回 -2）。

法线约定：n 由 A 指向 B。速度兜底取 A 相对 B 的来流方向 `A.v−B.v`，
冲量两侧按同一约定符号使用即自洽。

## 4. 双曲线定阈值（20 万接触，double 真值对照 float 观测）

`tools/sweep.c` 用 double 生成几何真值再 cast 成 float 喂被测算法，混合
普通尺度（L≤10）与大坐标消去（L 至 1e6）两类场景。沿对数阈值网格统计：

- **误判 misjudged**（蓝，上升）：真值 d_true ≥ 1e-4 几何可靠，float 却
  退化走兜底——阈值越大越多正常碰撞被误判；
- **坏法线 bad_normal**（橙，下降）：走了位置法线但与真值夹角 >30°，
  或 eps=0 时产出 NaN——阈值越小越多崩溃/错误法线。

关键数据（`results/sweep.csv`）：

| eps_t | 误判 | 坏法线 | 真崩溃 |
|---|---|---|---|
| 0（无守卫） | 274 | 820 | **9997** |
| 1e-7 | 0 | 4190 | 0 |
| **1e-4（选定）** | **274** | **820** | **0** |
| 2.15e-4（交叉网格点） | 7133 | 633 | 0 |
| 1e-1 | 61978 | 0 | 0 |

- 崩溃只在 eps=0 存活：任何正阈值都能消灭 0/0，但阈值必须高过噪声地板
  才能同时消灭"非零却错误"的法线。
- 交叉点在 1e-4～2e-4 之间。取 **eps_t = 1e-4**（`EPS_T_DEFAULT`）：
  位于噪声地板 1.2e-7 与交叉点之间偏保守一侧，相对交叉点用 274 次误判
  换掉约 600 次坏法线的余量；普通尺度下面板误判为 0——这个阈值只在
  大坐标消去场景才介入。
- 取舍是场景相关的：若仿真坐标长期在 L~1e6，应按 ULP(L) 上调阈值；
  L≤10 的常规场景 1e-4 已远离噪声地板。

## 5. 重合输入验收

球心完全重合的用例必须给出非零单位法线，自检逐条打印球编号与兜底方式：

```
[coincidence] minimal_coincident_opposite_velocity: balls=A,B mode=velocity-fallback normal=(1,0,0) non_zero_unit=1
[coincidence] counterexample_coincident_equal_velocity: balls=A,B mode=fixed-axis-fallback normal=(1,0,0) non_zero_unit=1
```

崩溃或输出零向量即 `all_ok=0`、进程非零退出。
