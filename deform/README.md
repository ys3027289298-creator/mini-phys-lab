# deform

柔性球接触：落地压扁 → 回弹。弹簧-阻尼接触模型，C 实现，`make test` 跑
单元自检 + 三档落高标定验证。

## 形变量表达式

接触期间沿法线方向的压缩量为 `x ≥ 0`（球竖向直径被压短的米数），
接触动力学：

```
m x'' = -k x - c x'        （x > 0 时；否则自由飞行）
```

- `k` 材料刚度 N/m，`c` 阻尼 N·s/m
- `wn = sqrt(k/m)`，阻尼比 `zeta = c / (2 m wn)`
- 最大压缩（最大压扁瞬间 `x' = 0`，能量守恒加阻尼修正）：

```
x_max = A(zeta) · v_impact / wn = A(zeta) · v_impact · sqrt(m / k)
A(zeta) = exp(-zeta · atan(sqrt(1-zeta²)/zeta) / sqrt(1-zeta²))
```

无阻尼极限 `A(0)=1`，退化为纯能量平衡 `½kx² = ½mv²`。

**最大压扁比**（视频逐帧量的量，竖向缩短量 ÷ 静止直径 `2r`）：

```
epsilon = x_max / (2r)
```

## 反解刚度

固定阻尼比，用目标压扁比反解：

```
k = m · (A(zeta) · v_impact / (2 r · epsilon_target))²
```

反弹由同一振子的恢复系数决定，与刚度解耦：

```
e = exp(-pi zeta / sqrt(1-zeta²))，   h_bounce = e² · h_drop
```

扫 `k` 时保持 `e` 固定，压扁可见度随 `1/sqrt(k)` 变化，反弹高度曲线基本
水平——这就是「果冻 ↔ 铅球」的取舍区间，见
`../artifacts/stiffness_tradeoff.png`（数据 `stiffness_sweep.csv`，
`make sweep` 重新生成）。软端反弹略降是接触时间变长、重力冲量变大的
真实效应，不是数值误差。

## 视频实测值与三档落高验证

逐帧量得的最大压扁比集中在 `src/video_measurements.h`，三档落高
0.5 / 1.0 / 2.0 m。程序用中档落高反解 `k`，再跑全部三档，每档打印
仿真压扁比、视频实测值、相对偏差；偏差超过 10% 判 FAIL 并额外打印
`actual_squash` 与 `video_squash` 两个数。

> **注意：参考视频不在工作区内。** 任务所述皮球视频未随工程提供，
> 机器上唯一的 `~/dock_take.mp4` 经核对是无关的 IDE 停靠录屏。
> `video_measurements.h` 当前是按 `sqrt(h)` 缩放规律构造的工程占位值，
> 文件顶部有醒目标注。真视频到位后的替换流程：
>
> 1. `ffmpeg -i ref.mp4 -vsync 0 frames/f%04d.png`
> 2. `python3 tools/measure_video.py frames --rgb R G B --out m.csv`
>    （颜色阈值用单帧先调；`python3 tools/smoke_measure.py` 是该工具
>    的合成球自测，已通过）
> 3. 把三次冲击的 epsilon 填回 `video_measurements.h`
> 4. `make test` 自动用新值重新标定、验证；超差会按要求报两个数

## 视觉形变

`deform_scales()` 给体积（面积）守恒的 squash & stretch 缩放：
竖向 `sy = 1 - epsilon`，横向 `sx = 1/sqrt(sy)`。

## 文件

| 文件 | 内容 |
|---|---|
| `src/deform.h` | 模型、解析式、反解刚度、时间步进仿真 |
| `src/deform.c` | 单元自检 + 标定 + 三档落高验证 + sweep CSV |
| `src/video_measurements.h` | 视频实测压扁比（**当前为占位值**） |
| `tools/measure_video.py` | 逐帧像素测量辅助（真视频到位后使用） |
| `tools/smoke_measure.py` | 测量工具的合成球烟雾测试 |
| `tools/plot_tradeoff.py` | 两条取舍曲线出图 |
