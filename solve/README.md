# solve

`make test` 跑原有单接触自检（输出末行 `all_ok=1`）。

## 30 球接触求解剖析基准

题述“30 球 / 8 ms 接触求解”的可计时重建与优化验证都在
[`src/bench_solve.c`](src/bench_solve.c)，完整报告见
[`../artifacts/contact_solve_report.md`](../artifacts/contact_solve_report.md)。

```bash
make bench          # 构建 bench_solve
make bench-test     # 基准健全性检查（接触数、有限性）
make bench-run      # 无参 = scenario + 三项分解 + 最终判定；达标退出码 0
./bench_solve help  # 所有模式
```

数据与图可复现：

```bash
taskset -c 2 python3 run_analysis.py   # 重测，写 ../artifacts/{*.csv,results.json,txt}
python3 make_curves.py                 # 重绘 ../artifacts/fig_*.png
```

结论（本机 Ryzen 9700X，绑核中位）：现状固定 4096 轮实测 8.43 ms，
其中接触枚举 0.03%、每轮重复派生计算 45.1%、冲量核心 54.9%；
T(K) ≈ 48.5 + 2.06K µs（R²=0.995）。改用“派生数据只算一次 +
相邻轮最大法向冲量变化 < ε 即停”，取 ε=0.01（由位置残差 1e-3 反解并经
阈值扫描整定），冷重载帧 5 轮收敛，求解中位 8.7 µs / p95 19.1 µs，
2400 帧最坏位置残差 2.24e-4 < 1e-3，`overall=1`。
