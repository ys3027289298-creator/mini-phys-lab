# mini-phys-lab

刚体碰撞与响应基础模块，C 实现，每目录独立构建可自测。

| 目录 | 内容 |
|---|---|
| `collide/` | 球与盒的接触生成、法线、穿透深度 |
| `solve/` | 法向冲量与切向摩擦、恢复系数 |
| `integrate/` | 半隐式欧拉、位置修正、不穿地 |
| `angular/` | 碰撞总角动量预算、千分之一守恒校验、力臂精度反解与耗时基准 |

角动量不等式推导、力臂误差曲线与精度取舍见
[`artifacts/angular_momentum_report.md`](artifacts/angular_momentum_report.md)。
