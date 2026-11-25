
/* sweep: 阈值 eps_t 扫描台，用双曲线定阈值。
 *
 * 每个样本用 double 生成几何真值，再把位置 cast 成 float 喂给被测算法，
 * 因此 float 减法消去（大坐标处）会真实发生。
 *
 * 对每个 eps 统计接触样本上的（几何可靠群体锚定为 d_true >= EPS_REF）：
 *   misjudged   真值可定法线(d_true>=EPS_REF)，float 却退化走兜底 —— 随 eps 上升
 *   bad_normal  走位置法线，但与真值夹角 >30deg，或产出 NaN（eps=0 的 0/0）—— 随 eps 下降
 *   crashes     d_flt==0 且阈值未覆盖（守卫 eps<=0 时才发生的真除零）
 *
 * 场景混合：90% 普通尺度(L<=10)，10% 大坐标平移(L 至 1e6)制造消去噪声；
 * 5% 精确重合；10% 两球同速（触发固定轴路径，不参与法线误差统计）。
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include "collide_ss.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define N 200000
#define N_GRID 26          /* 0 + 25 个对数点 */
#define ANGLE_BAD (30.0 * M_PI / 180.0)
#define EPS_REF 1e-4       /* 几何可靠群体锚点：d_true 高于此值即视为本可定法线 */

typedef struct {
    double dx, dy, dz;     /* 真值相对位置（A->B） */
    float  fx, fy, fz;     /* cast 成 float 后的相对位置 */
    double d_true;
    float  d_flt;
    int    exact_coincident;
    int    same_velocity;
    int    large_scale;     /* 大坐标平移诱发的消去噪声样本 */
} Pair;

static uint64_t rng = 0x123456789abcdef0ULL;
static double urand(void) {            /* (0,1) */
    rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
    return (rng >> 11) * (1.0 / 9007199254740992.0);
}
static double grand(void) {            /* Box-Muller 标准正态 */
    double u1 = urand(), u2 = urand();
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

int main(int argc, char** argv) {
    const char* out = (argc > 1) ? argv[1] : "results/sweep.csv";
    Pair* ps = malloc(sizeof(Pair) * N);
    int contacts = 0;

    for (int i = 0; i < N; i++) {
        Pair* p = &ps[i];
        /* 球心距：对数均匀，1e-9 .. 1.5，覆盖深穿透到分离 */
        double logd = -9.0 + urand() * (log10(1.5) + 9.0);
        double d = pow(10.0, logd);
        if (urand() < 0.05) d = 0.0;                /* 5% 精确重合 */

        double gx = grand(), gy = grand(), gz = grand();
        double gl = sqrt(gx * gx + gy * gy + gz * gz);
        p->exact_coincident = (d == 0.0);
        p->d_true = d;
        p->dx = d * gx / gl; p->dy = d * gy / gl; p->dz = d * gz / gl;

        /* 10% 概率整体平移到大坐标 L，诱发 float 减法消去 */
        double shift = 0.0;
        p->large_scale = 0;
        if (urand() < 0.10) {
            double L = pow(10.0, urand() * 6.0);   /* 1 .. 1e6 */
            shift = L;
            p->large_scale = 1;
        }
        float ax = (float)(shift), ay = (float)(shift * 0.5), az = 0.0f;
        float fx = (float)(shift + p->dx) - ax;
        float fy = (float)(shift * 0.5 + p->dy) - ay;
        float fz = (float)(p->dz) - az;
        p->fx = fx; p->fy = fy; p->fz = fz;
        p->d_flt = sqrtf(fx * fx + fy * fy + fz * fz);
        p->same_velocity = (urand() < 0.10);

        if (p->d_flt < 1.0f) contacts++;            /* rsum=1，接触事件 */
    }

    double grid[N_GRID];
    grid[0] = 0.0;
    for (int k = 1; k < N_GRID; k++)
        grid[k] = pow(10.0, -9.0 + (k - 1) * (8.0 / (N_GRID - 2))); /* 1e-9..1e-1 */

    FILE* f = fopen(out, "w");
    fprintf(f, "eps,contacts,misjudged,bad_normal,crashes,fallbacks,"
               "misjudged_local,bad_normal_local,misjudged_large,bad_normal_large\n");

    for (int k = 0; k < N_GRID; k++) {
        double eps = grid[k];
        long misjudged = 0, bad = 0, crashes = 0, fallbacks = 0;
        long mis_loc = 0, bad_loc = 0, mis_lrg = 0, bad_lrg = 0;

        for (int i = 0; i < N; i++) {
            Pair* p = &ps[i];
            if (p->d_flt >= 1.0f) continue;
            Ball A = { {0, 0, 0}, {1.0f, 0.0f, 0.0f}, 0.5f };
            Ball B = { {p->fx, p->fy, p->fz},
                       p->same_velocity ? A.v : (V3){ -1.0f, 0.0f, 0.0f }, 0.5f };
            Contact c;
            balls_contact(&A, &B, (float)eps, &c);

            int guarded = (p->d_flt < eps);
            if (!guarded && p->d_flt == 0.0f) crashes++;   /* 真 0/0 未被覆盖 */
            if (c.mode != MODE_POSITION) {
                fallbacks++;
                /* 误判：几何上可靠（真值高于锚点、非精确重合），却被阈值兜走 */
                if (p->d_true >= EPS_REF) {
                    misjudged++;
                    if (p->large_scale) mis_lrg++; else mis_loc++;
                }
            } else if (!p->exact_coincident && p->d_true > 0.0) {
                double ntx = p->dx / p->d_true, nty = p->dy / p->d_true,
                       ntz = p->dz / p->d_true;
                double dot = ntx * c.normal.x + nty * c.normal.y + ntz * c.normal.z;
                if (!isfinite(dot)) { bad++; if (p->large_scale) bad_lrg++; else bad_loc++; continue; }
                if (dot < -1.0) dot = -1.0;
                if (dot > 1.0) dot = 1.0;
                if (acos(dot) > ANGLE_BAD) {
                    bad++;
                    if (p->large_scale) bad_lrg++; else bad_loc++;
                }
            }
        }
        fprintf(f, "%.3g,%d,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld\n",
                eps, contacts, misjudged, bad, crashes, fallbacks,
                mis_loc, bad_loc, mis_lrg, bad_lrg);
    }
    fclose(f);
    free(ps);
    printf("wrote %s (contacts=%d / %d pairs)\n", out, contacts, N);
    return 0;
}
