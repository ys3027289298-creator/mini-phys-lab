# angular

碰撞总角动量预算：三维刚体法向冲量（含转动有效质量与 `r×J` 角冲量）、
碰撞前后总角动量独立校验（千分之一红线）、力臂精度扫描与耗时基准。

符号约定、不等式推导、实测曲线与精度决策见
[`../artifacts/angular_momentum_report.md`](../artifacts/angular_momentum_report.md)。

```sh
make test       # 自检：手算值 + 违规上报 + 2 万随机偏心碰撞
make bench      # 构建基准 am_bench
make bench-run  # 跑基准并把 am_*.csv 写入 ../artifacts
make plot       # 由 CSV 重绘 ../artifacts/am_*.png
```

违规时程序输出一行，含角动量变化量与当时的偏心距：

```text
AM_VIOLATION case=... dL_kgm2/s=... rel=... limit=1e-03 ecc_max_m=... eccA_m=... eccB_m=...
```

关键结论（SI，典型偏心算例 J≈5.97 N·s、L_ref=3.0 kg·m²/s）：
力臂误差每 1 mm 增加约 1.99‰ 相对角动量误差，反解最大力臂误差 ≈0.50 mm；
随机包络下保守量化指标 ≈0.05 mm。double 全叉积实测偏差 ~1e-16、约 35 ns/次，
加逐次校验约 70 ns/次，故选「double 全叉积 + 常开校验」。
