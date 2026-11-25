
/* collide_ss: 球-球接触生成（normal + depth），含退化法线三级兜底。
 *
 * 法线约定：n 由 A 指向 B（沿 B.c - A.c）。
 *
 * 退化判据：球心距 d < eps_t 时，位置差无法确定法线方向。
 *   - 几何含义（等半径 r=0.5 时 rsum=1）：此时穿透深度 δ = 1-d > 1-eps_t，
 *     即 eps_t 是“仍允许用位置定法线”的最小球心距，等价地由预期最大
 *     重叠深度反解：eps_t = rsum - δ_max。
 *   - 数值含义：float32 在坐标尺度 L 处的分辨力 ULP(L) ~= 1.2e-7 * L，
 *     d 低于该尺度时减法消去，算出的 d 只是 0 或噪声，eps_t 必须高于
 *     这个噪声地板。
 *
 * 兜底顺序（每级都可能单独退化，故不能只做一级）：
 *   1. 位置差归一化（几何法线）
 *   2. 相对速度 vA-vB 归一化（运动学法线；两球同速时同样退化为零向量）
 *   3. 固定轴 (1,0,0)（永不退化，球心完全重合也给出非零单位法线）
 */
#ifndef COLLIDE_SS_H
#define COLLIDE_SS_H

typedef struct { float x, y, z; } V3;
typedef struct { V3 c; V3 v; float r; } Ball;
typedef struct { V3 normal; float depth; int hit; int mode; } Contact;

/* 兜底方式标记，供调用方/自检核对 */
enum { MODE_POSITION = 0, MODE_VELOCITY = 1, MODE_FIXED_AXIS = 2 };

/* 缺省阈值：高于 L~1 的 float 噪声地板(~1.2e-7)，对应 δ_max≈0.9999 */
#define EPS_T_DEFAULT 1e-4f

/* 球-球接触。eps_t<=0 时使用缺省阈值。
 * 返回 1=有接触, 0=分离, -1=非法半径, -2=内部不变量被破坏(未输出法线) */
int balls_contact(const Ball* a, const Ball* b, float eps_t, Contact* o);

#endif
