
/* ledger: 落球-弹起能量账（N/T/C/R 四项 + 残差反解 + 四档落高闭合自检）。
 *
 * 坐标 y 为球心高度，地面 y=0，接触位置 y=radius；势能零点取接触高度，
 * 故 E = m g (y-radius) + 1/2 m (vx^2+vy^2)。
 *
 * 四项一律按“带符号的能量变化贡献”记账（损失为负、注入为正）：
 *   N = -1/2 m [(vn-)^2-(vn+)^2] = -1/2 m (vn-)^2 (1-e^2)
 *   T = -1/2 m [(vt-)^2-(vt+)^2],  vt+ = max(|vt-|-mu|vn-|,0)
 *   C = + mg dy_corr + 1/2 m (|v_clamp+|^2-|v_clamp-|^2)
 *   R = dE - N - T - C                      （总能量变化反解）
 * 恒等式 N+T+C+R == dE 为构造保证；实质性闭合看 |R|/E0。
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MASS   1.0f
#define GRAV   9.81f
#define RADIUS 0.5f

typedef struct { float x, y, vx, vy; } Ball;
typedef struct { float e, mu, dt; } SimCfg;

typedef struct {
    int contact;
    float pen, dy;
    float vn_pre, vt_pre, vn_post, vt_post;
    float dE, N, T, C, R;
} Frame;

static float mech_E(const Ball* b, int include_tangent) {
    float tang = include_tangent ? b->vx * b->vx : 0.0f;
    return MASS * GRAV * (b->y - RADIUS)
         + 0.5f * MASS * (tang + b->vy * b->vy);
}

/* 与 solve.c 同一冲量模型、与 integrate.c 同一更新/修正顺序，并逐项留账。 */
static void step_ledger(Ball* b, const SimCfg* cf, Frame* f) {
    memset(f, 0, sizeof(*f));
    float E_start = mech_E(b, 1);

    /* 阶段 1：半隐式欧拉（先更新速度，再用新速度更新位置） */
    b->vy -= GRAV * cf->dt;
    b->y  += b->vy * cf->dt;

    /* 阶段 2：法向恢复 + 切向摩擦冲量 */
    if (b->y < RADIUS) {
        f->contact = 1;
        float vn = b->vy, vt = b->vx;
        f->vn_pre = vn;
        f->vt_pre = vt;
        if (vn < 0.0f) {
            float vn_post = -cf->e * vn;
            float u = fabsf(vt);
            float reduce = cf->mu * fabsf(vn);
            if (reduce > u) reduce = u;          /* 摩擦不能反向 */
            float scale = (u > 1e-9f) ? (u - reduce) / u : 0.0f;
            float vx_post = scale * vt;

            f->N = -0.5f * MASS * (vn * vn - vn_post * vn_post);
            f->T = -0.5f * MASS * (vt * vt - vx_post * vx_post);

            b->vy = vn_post;
            b->vx = vx_post;
            f->vn_post = vn_post;
            f->vt_post = vx_post;
        }

        /* 阶段 3：位置投影修正 + 穿透帧速度手术（在 pen>0 块内） */
        f->pen = RADIUS - b->y;
        f->dy = f->pen;
        b->y += f->pen;
        float vy_pre_clamp = b->vy;
        if (b->vy < 0.0f) b->vy = 0.0f;
        float dK_clamp = 0.5f * MASS * (b->vy * b->vy - vy_pre_clamp * vy_pre_clamp);
        f->C = MASS * GRAV * f->dy + dK_clamp;
    }

    float E_end = mech_E(b, 1);
    f->dE = E_end - E_start;
    f->R = f->dE - f->N - f->T - f->C;   /* 残差反解 */
}

typedef struct {
    float E0, E1, dE_sum, N, T, C, R;
    int frames, contact_frame;
} Account;

/* 从 H（接触高度以上的落高）静止释放，记到第一次弹起后的最高点。 */
static void run_drop(float H, float e, float mu, float dt, Account* a) {
    Ball b = { 0.0f, RADIUS + H, 0.0f, 0.0f };
    SimCfg cf = { e, mu, dt };
    memset(a, 0, sizeof(*a));
    a->E0 = mech_E(&b, 1);
    a->contact_frame = -1;
    int seen = 0, k;
    for (k = 0; k < 200000; k++) {
        Frame f;
        step_ledger(&b, &cf, &f);
        a->frames++;
        a->dE_sum += f.dE;
        a->N += f.N; a->T += f.T; a->C += f.C; a->R += f.R;
        if (f.contact && !seen) { seen = 1; a->contact_frame = k; }
        if (seen && k > a->contact_frame && b.vy <= 0.0f) break;
    }
    a->E1 = mech_E(&b, 1);
}

/* 差额最大项挑选：先 N/T/C 账实核对，三者都相符才挑 R。返回 0N 1T 2C 3R。 */
static int pick_worst(float dN, float dT, float dC, float R, float E0) {
    float tol = 1e-6f * E0;   /* float 版精确复核容差放宽到 1e-6 E0 */
    float mx = dN; int pick = 0;
    if (dT > mx) { mx = dT; pick = 1; }
    if (dC > mx) { mx = dC; pick = 2; }
    if (mx > tol) return pick;
    (void)R;
    return 3;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* ---- 公式单点核对：vn=-10,e=.5 -> Ln=37.5 J ---- */
    float vn = -10.0f, e = 0.5f;
    float Ln = 0.5f * MASS * vn * vn * (1.0f - e * e);
    check("normal_formula_37_5", fabsf(Ln - 37.5f) < 1e-4f);

    /* vt=4, |vn|=1, mu=1 -> reduce=1, vt'=3 -> Lt=.5(16-9)=3.5 J */
    float vt = 4.0f, mu = 1.0f;
    float reduce = fminf(mu * 1.0f, vt);
    float vt2 = vt - reduce;
    float Lt = 0.5f * MASS * (vt * vt - vt2 * vt2);
    check("tangent_formula_3_5", fabsf(Lt - 3.5f) < 1e-4f);

    /* ---- “R≡0 且 T≡0”诊断：solve.c 自检同款 v=(4,-1), e=0, mu=1 ---- */
    {
        float nvn = -1.0f, nvt = 4.0f, ne = 0.0f, nmu = 1.0f;
        float vn2 = -ne * nvn;
        float red = fminf(nmu * fabsf(nvn), fabsf(nvt));
        float tv2 = fmaxf(fabsf(nvt) - red, 0.0f);
        float lN = 0.5f * (nvn * nvn - vn2 * vn2);       /* 0.5 */
        float lT = 0.5f * (nvt * nvt - tv2 * tv2);       /* 3.5 */
        float dE_full = -(lN + lT);                      /* -4.0 */
        float dE_vert = 0.5f * vn2 * vn2 - 0.5f * nvn * nvn; /* -0.5 */
        float RA = dE_full - (-lN) - (-lT);               /* 完整账: 0 */
        float RB = dE_full - (-lN);                       /* 全能量漏 T: -3.5 */
        float RC = dE_vert - (-lN);                       /* 截断竖直口径: 0, T=0 */
        check("demo_A_R_zero", fabsf(RA) < 1e-6f);
        check("demo_B_R_absorbs_tangent", fabsf(RB + 3.5f) < 1e-5f);
        check("demo_C_R_zero_and_T_zero", fabsf(RC) < 1e-6f);
        check("demo_tangent_impulse_real", fabsf(tv2 - 3.0f) < 1e-6f);
        printf("demo_three_ledgers: A(N=%.2f,T=%.2f,R=%.2e) "
               "B(N=%.2f,T=0,R=%.2f) C(N=%.2f,T=0,R=%.2e)\n",
               -lN, -lT, RA, -lN, RB, -lN, RC);
    }

    /* ---- picker 规则自检 ---- */
    check("picker_picks_normal",    pick_worst(24.8f, 0.0f, 0.0f, 25.9f, 47.0f) == 0);
    check("picker_picks_tangent",   pick_worst(0.0f, 7.1f, 0.0f, 8.2f, 47.0f) == 1);
    check("picker_picks_corr",      pick_worst(0.0f, 0.0f, 0.46f, 0.64f, 47.0f) == 2);
    check("picker_clean_picks_R",   pick_worst(1e-9f, 1e-9f, 1e-9f, 0.03f, 9.8f) == 3);

    /* ---- 四档落高能量账闭合 ---- */
    const float Hs[4] = { 1.0f, 2.0f, 4.0f, 8.0f };
    int all_close = 1;
    puts("H,frames,E0,E1,N,T,C,R,resid_pctE0,identity_pctE0,pass");
    for (int i = 0; i < 4; i++) {
        Account a;
        run_drop(Hs[i], 0.5f, 0.2f, 1e-3f, &a);
        float resid_pct = 100.0f * fabsf(a.R) / a.E0;
        float sum4 = a.N + a.T + a.C + a.R;
        float ident_pct = 100.0f * fabsf(sum4 - (a.E1 - a.E0)) / a.E0;
        int pass = (resid_pct < 1.0f) && (ident_pct < 1.0f);
        if (!pass) all_close = 0;
        printf("H=%.0f frames=%d E0=%.3f E1=%.3f N=%.4f T=%.4f C=%.5f R=%.5f "
               "resid_pct=%.4f identity_pct=%.2e pass=%d\n",
               Hs[i], a.frames, a.E0, a.E1, a.N, a.T, a.C, a.R,
               resid_pct, ident_pct, pass);
        /* 弹起能量应在理想 e^2=0.25 附近（修正注入 + 积分漂移，带宽 0.18..0.32） */
        check("height_bounce_ratio_sane",
              a.E1 / a.E0 > 0.18f && a.E1 / a.E0 < 0.32f);
    }
    check("four_heights_close_1pct", all_close);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
