# 碰撞系统故障注入研究（fi/）

对 `collide/`、`solve/`、`integrate/` 三个模块做主动故障注入，验证碰撞系统的健壮性。
本目录不改动被测源码，通过"包含源码 + `#undef main`"的适配器直接调用其 `static` 内部函数。

复现：`./fi/run_all.sh`（依赖 gcc / gcov / python3 / scipy；退出码 0 表示 9 用例全通过）。

---

## 1. 三类故障与判据表达式

每个判据是若干**可计算的比较式**的合取，检查项按"浅 → 深"排列。这个顺序是有意的：
在**有故障的调用上仍然通过的浅层检查**，就是"被掩盖的那一层"。

记法：`r` 返回码；`O` 纯输出参数（`Contact c`、`float loss`）；`S` 权威状态（`V3 v`、`Body b`）；
`fin(x)` 非 NaN/±Inf；`deg` 注入的退化数值。

```
C_num   = in_codes(r)          # r ∈ 模块文档化返回码
        ∧ finite_output(O)     # 输出数值有限、范围合法
        ∧ reject_deg           # r < 0  当且仅当注入的是退化数

C_state = defined(O/S)         # 所有上报字段有限、非 NaN
        ∧ cleared(O)           # miss/reject 时槽位不残留上一帧的值
        ∧ atomic(S)            # reject 时权威状态逐位不变

C_time  = ret_eq               # 两序返回码一致
        ∧ fin_eq               # 两序输出均有限
        ∧ val_eq               # 每个可观测量逐位相等（含 maxdiff == 0）
```

**用例判定**：本类判据通过 = PASS。某故障使自己那类判据通过、却在同一次调用上被
另一类判据抓到 = **被掩盖的故障**，记录"哪一层失效"。

### 为什么 `cleared` 是语义必需项，而非苛求

任务判定式 `上帧残留接触...应被清除` 直接要求它。更根本地：`Contact`/`loss` 是**纯输出参数**，
调用方在 `r=0`（miss）或 `r<0`（reject）之后照样读取它们。若不清零，调用方拿到的就是上一帧的
穿透深度与法线——这正是状态污染。三模块里 `solve.resolve` 与 `integrate.step` 都在拒绝路径上
正确地把输出清零/保持原子；唯独 `collide.test` 没有，因此这是真实的模块间不一致，而非判据造出来的。

### 注入点选择（不是随便选）

每个注入点都落在源码中一条**真实判定路径**上，并被 gcov 证实能到达目标代码：

| 用例 | 注入 | 命中源码位置 | 该点覆盖的路径（节选） |
|---|---|---|---|
| N1 | 球半径 = 0 | `collide.c:29` `if (s->r <= 0.0f) return -1` | `entry:test`, `guard:radius` |
| N2 | 接触法线 = (0,0,0) | `solve.c:34` `scale = (tl>1e-9f)? ... : 0` | 13 条 resolve 路径 |
| N3 | dt = +Inf | `integrate.c:10` `if (dt < 0.0f …)` — Inf 穿透该守卫 | 10 条 step 路径 |
| S1 | 上帧 Contact 种子 | `collide.c:34-35` 只重置 `hit/depth` | 10 条路径（均不含 `decision:surface`） |
| S2 | 上帧 loss 种子 | `solve.c:16` 提前 return，不写 `*out_normal_loss` | 仅 `entry:resolve`,`guard:restitution` |
| S3 | 上帧 prev-position 种子 | `integrate.c:10` 提前 return，`px/py` 不动 | 8 条路径 |
| T1 | 4 路调用置换 | `collide.test` 全部路径 | 26/26 collide 路径 |
| T2 | 3 路调用置换 | `solve.resolve` 全部路径 | 13/13 solve 路径 |
| T3 | 3 路调用置换 | `integrate.step` 全部路径 | 10/10 integrate 路径 |

---

## 2. 覆盖反解注入点数

**覆盖单元（"路径"的严格代理）** = 可执行源码行 ∪ 分支臂，按函数归属。
范围为三个模块的功能函数：`clampf/closest/test`、`resolve`、`step`。
**排除**模块自带自测（`check` + 内嵌 `main`）：没有任何注入点驱动它们，计进去只会无差别压低所有点的覆盖率。

要诚实说明两点：
1. 这是**可达单元覆盖**，不是"真路径"（分支交错组合）枚举；
2. 命名字段的"路径"粒度（51 个语义区域）太粗——3 个点就达 96%，无法区分优劣。
   因此选点用**单元粒度**（124 个单元），路径粒度只用于向人解释"覆盖了哪些路径"。

**候选点 → 覆盖集合**（`results/coverage_matrix.json`、`point_paths.json`）。
**时间曲线**：每个注入点是一个独立进程（这正是套件的运行方式），实测单进程 ≈ 0.76 ms，
近乎固定开销，故耗时基本随 k 线性（中位进程耗时见 `point_times.json`，另附进程内单次调用耗时）。

**覆盖率曲线**（贪心按绝对增量，`results/tradeoff.txt`）：

```
k= 1  collide/T1         +52  52/124  41.9%
k= 2  solve/T2           +26  78/124  62.9%
k= 3  integrate/T3       +17  95/124  76.6%   <- 拐点/膝点
k= 4  collide/P_insideYl  +4  99/124  79.8%
...
k= 9  (强制 3/类)            99/124  79.8%
k=11  (精确最优)            112/124  90.3%   <- 90% 最小点
k=17                         118/124  95.2%
k=23                         124/124 100.0%
```

**两个曲线定点**：
- **拐点** k=3（76.6%）：三个 T 用例各覆盖一个完整模块，边际收益 17~52 单元/点；
  此后每点只 +1~4 单元，单位时间收益从 ~25 掉到 ~1.4 单元/ms。**投入产出比在此崩塌**。
- **90% 精确最小点** k=11（用 ILP 求，非贪心上界）：
  `results/optimal_sets.json`。
- 100% 需 k=23；95% 需 k=17。

**取舍结论**：
- 若目标是**发现缺陷**（本任务），k=9 的三类各三点已经覆盖了 90% 路径所不能区分的**判定语义**——
  9 个用例全部落在不同的语义决策上，且暴露了 4 个真实缺陷（见 §3）。此时再加点只增耗时、不增语义。
- 若目标是**行/分支覆盖数字**，则强制 3/类只有 79.8%，需 **加 5 点、k=14 才到 90.3%**；
  纯覆盖率最优则是 **k=11**（可少于 9+5，因为它不要求每类都有点）。
- 因此：**9 个用例保留为语义核心；要补覆盖率就按 `optimal_sets.json` 补 5 个**，
  而不是按类别均摊地堆点数。

---

## 3. 九个用例的通过/失败结论

| 用例 | 类 | 注入 | 结论 | 首个失败检查 | 失效的上层判据层 |
|---|---|---|---|---|---|
| N1 | 数值 | 半径 0 | **PASS** | — | 见下方"被掩盖的故障" |
| S1 | 状态 | 上帧 Contact | **FAIL** | `cleared` | `defined` 层失明 |
| T1 | 时序 | 4 路置换 | **PASS** | — | — |
| N2 | 数值 | 法线 (0,0,0) | **FAIL** | `reject_deg` | `in_codes`、`finite_output` 层失明 |
| S2 | 状态 | 上帧 loss | **FAIL** | `cleared` | `defined` 层失明 |
| T2 | 时序 | 3 路置换 | **PASS** | — | — |
| N3 | 数值 | dt = +Inf | **FAIL** | `finite_output` | `in_codes` 层失明 |
| S3 | 状态 | 上帧 prev-pos | **PASS** | — | — |
| T3 | 时序 | 3 路置换 | **PASS** | — | — |

**PASS 5 / FAIL 4**。各失败的具体证据：

- **S1**：`r=0` 的 miss 返回后 `hit=0,depth=0`，但 `normal=(0.7,-0.7,0.1)`、`point=(7,8,9)`
  仍是上一帧的值——`collide.test` 只在成功路径写全 `Contact`。
- **S2**：reject（`restitution=1.5`）返回 `r=-1`，`loss` 仍为注入前种子 `0.123`。
  `resolve` 在所有 reject 路径上都不写 `*out_normal_loss`。
- **N2**：零长度法线使 `vn=0`（不满足 `vn>0` 的 separating 分支），`resolve` 返回 **1（成功）**，
  但速度 `v` 纹丝不动、`loss=0`——**静默 no-op 被当作成功**。`in_codes`（r=1 合法）与
  `finite_output`（输出全有限）两层都放行了，直到 `reject_deg` 才拦下。
- **N3**：`dt=+Inf` 满足 `dt<0` 为假而穿透守卫；`y/=NaN`、`vx=NaN`、`vy` 被地板修正成 0，
  得到"部分 NaN"的混合状态。`in_codes`（r=1）先放行，`finite_output` 才拦下。

### 被掩盖的故障（mission 明确要求指出）

1. **N1 数值故障被自身判据放过**：`C_num` 全通过（r=-1 合法、无 NaN），但同一调用上
   `C_state.cleared` 失败——reject 后 `Contact` 仍带着上一帧的 `normal/point`。
   即**数值退化故障本身被"数值判据"完全遮住，是"状态判据"抓到的**。
   这也解释了 S1 的根因：`collide.test` 的任何非成功返回都不定义输出槽位。
2. **N2 的失效层**：`in_codes` 与 `finite_output` 两层都在"有故障的调用"上通过，
   故障只在最深的 `reject_deg` 才暴露——**浅层守卫对"语义无效但数值合法"的输入是盲的**。
3. **N3 的失效层**：`in_codes` 放行（r=1 在文档码内）；`finite_output` 才拦下。

### 时序类为何三个全 PASS（以及它掩盖了什么）

`T1/T2/T3` 用的是**相互独立的输出槽位**。三个模块的 `test/resolve/step` 对独立槽位是纯函数、
无全局可变状态，因此任何置换都逐位相等——**PASS 是真实结论**。
但这恰恰掩盖了时序脆弱性：一旦**多个相处共享同一输出槽位**，`collide.test` 的残留就变成顺序依赖。
本目录提供补充诊断 `collide X_alias` 复现之：

```
共享槽 [hit, miss]  -> hit=0 depth=0
共享槽 [miss, hit]  -> hit=1 depth=0.1      equal=0
```

即**同一组调用只因顺序不同而结果不同**。这正是"时序错乱被状态污染掩盖"的实例：
时序判据在独立槽位下无法暴露它，必须靠状态判据（S1 的 `cleared`）从根上消除残留。
换句话说，**修掉 `cleared` 缺陷，时序类的隐患连带消失**——这也是为什么状态判据是本套件的核心。

---

## 4. 目录与产物

```
fi/
  adapters/{collide,solve,integrate}_adapter.c   注入适配器（含源码、bench、list）
  build.sh                                       构建普通版 + gcov 版
  run_all.sh                                     一键复现
  tools/covtool.py   逐点覆盖率采集（GCOV_PREFIX 隔离 + gcno 镜像）、声明集、计时
  tools/criteria.py  三类判据与分层判定
  tools/paths.py     源码行 ↔ 命名路径词典
  tools/run_suite.py 跑 9 用例、出结论与掩盖分析
  tools/tradeoff.py  覆盖/时间曲线、精确最小 k（ILP）、定点决策
  results/verdicts.json         9 用例完整判定
  results/coverage_matrix.json  每个候选点的覆盖单元集合
  results/point_paths.json      每个候选点的命名路径
  results/coverage_curve.json / tradeoff.txt   两条曲线与报告
  results/optimal_sets.json     k=90/95/100 的最优点集
  results/decision.json         最终决策参数
```

## 5. 局限

- 覆盖单元是**行 + 分支臂**，不是分支交错组合的真路径；"90% 路径"按此口径。
- 计时是单机单次进程级中位数，固定 exec 开销主导；进程内单次调用耗时另表列出，二者都给出。
- 候选点池（52 个）由人工按源码判定路径构造，已用"声明集 ⊆ 并集"断言证明池子能 100% 覆盖，但不保证是**最小**候选池。
- 命名路径词典按本 commit 行号固定，源码行变动需同步更新（`paths.py`）。
