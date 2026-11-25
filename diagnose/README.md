# diagnose

`make test` 运行日志归因与周期公式的自检。

分析真实日志：

```sh
python3 -I audit_log.py contact.log --fps 60 --speed 3.0 --radius 0.5
```

单行结构化字段越完整，能排除的路径越多；最少应包含 `frame ball depth normal radius center vn arm awake solver_called why j`。

## 日志粒度取舍

设一帧物理耗时为 `T`，写一条日志的耗时为 `c_g`，真正落盘的帧比例为 `p_g`，则：

```text
FPS_g / FPS_0 = T / (T + p_g c_g)
帧率损失 L_g  = p_g c_g / (T + p_g c_g)
```

周期异常约每 `N=30~60` 帧一次时，异常触发日志的 `p_g≈1/N`。若 `ρ=c_g/T=0.20`，每帧全量日志损失 `ρ/(1+ρ)=16.7%`；只在异常帧写则损失 `(ρ/N)/(1+ρ/N)`，约 `0.33%~0.66%`。

推荐粒度：

1. 每帧只写固定大小的内存 ring（保留最近 32~64 帧，无文件 I/O）。
2. 仅在 `contact && (j==0 || |Δw|==0)` 时输出一行紧凑归因记录，并转储 ring。
3. 每 30/60 帧输出一条 heartbeat。

四个独立判定位 `geometry_ok / vn_gate / arm_gate / sleep_gate` 可区分 4 类单因；位组合最多给出 `2^4-1=15` 种非空多因签名。若固定等间隔采样而不是异常触发，采样间隔 `g` 必须满足 `g < N_min/2`；对“几十帧”的周期，`g≤10` 才不至于混叠。
