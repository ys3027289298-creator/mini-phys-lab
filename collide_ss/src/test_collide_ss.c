
/* test_collide_ss: 球-球退化法线自检。
 * 打开 FE_DIVBYZERO/FE_INVALID：任何残留的 0/0 立即 SIGFPE，而非静默 NaN。
 *
 * 球心完全重合的用例逐条打印球编号与所用兜底方式；
 * 任一用例崩溃或输出零向量，final verdict 为 FAIL。 */
#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <fenv.h>
#include <string.h>
#include "collide_ss.h"

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static int is_unit_nonzero(const Contact* c) {
    float l2 = c->normal.x * c->normal.x + c->normal.y * c->normal.y
             + c->normal.z * c->normal.z;
    return l2 > 0.0f && isfinite(l2) && fabsf(sqrtf(l2) - 1.0f) < 1e-5f;
}

static const char* mode_name(int mode) {
    switch (mode) {
        case MODE_POSITION:   return "position";
        case MODE_VELOCITY:   return "velocity-fallback";
        case MODE_FIXED_AXIS: return "fixed-axis-fallback";
        default:              return "unknown";
    }
}

/* 球心完全重合的专项验收：必须给出非零单位法线，并记录球编号与兜底方式。
 * 接触写回调用方的 Contact，便于后续断言具体模式。 */
static void coincidence_acceptance(const char* tag, const Ball* A, const Ball* B,
                                   Contact* c) {
    int rc = balls_contact(A, B, EPS_T_DEFAULT, c);
    int ok = (rc == 1) && is_unit_nonzero(c);
    printf("[coincidence] %s: rc=%d balls=A,B mode=%s normal=(%g,%g,%g) "
           "depth=%g non_zero_unit=%d\n",
           tag, rc, mode_name(c->mode), c->normal.x, c->normal.y, c->normal.z,
           c->depth, ok);
    check(tag, ok);
}

int main(void) {
    feenableexcept(FE_DIVBYZERO | FE_INVALID);
    Contact c;

    /* 原始输入（原点附近，d=0.02，可分辨）：走几何法线，不崩 */
    Ball A = { {0.0f, 0.0f, 0.0f}, { 1.0f, 0.0f, 0.0f}, 0.5f };
    Ball B = { {0.02f, 0.0f, 0.0f}, {-1.0f, 0.0f, 0.0f}, 0.5f };
    check("original_near_origin", balls_contact(&A, &B, EPS_T_DEFAULT, &c) == 1);
    check("original_position_mode", c.mode == MODE_POSITION);
    check("original_normal_x", fabsf(c.normal.x - 1.0f) < 1e-5f);
    check("original_depth", fabsf(c.depth - 0.98f) < 1e-5f);
    check("original_unit", is_unit_nonzero(&c));

    /* 同一输入搬到大坐标：float 相减消去，d_flt=0，相对速度仍可定法线 */
    Ball A2 = { {1.0e6f, 0.0f, 0.0f}, { 1.0f, 0.0f, 0.0f}, 0.5f };
    Ball B2 = { {1.0e6f + 0.02f, 0.0f, 0.0f}, {-1.0f, 0.0f, 0.0f}, 0.5f };
    check("large_offset_no_crash", balls_contact(&A2, &B2, EPS_T_DEFAULT, &c) == 1);
    check("large_offset_velocity_mode", c.mode == MODE_VELOCITY);
    check("large_offset_unit", is_unit_nonzero(&c));
    /* 来流 A 相对 B 为 (2,0,0)，约定 n 由 A 指向 B，故为 +x */
    check("large_offset_normal_x", fabsf(c.normal.x - 1.0f) < 1e-5f);

    /* 最小输入：球心完全重合、等半径、速度反向 -> 速度兜底 */
    Ball C = { {3.0f, 4.0f, 5.0f}, { 1.0f, 0.0f, 0.0f}, 0.5f };
    Ball D = { {3.0f, 4.0f, 5.0f}, {-1.0f, 0.0f, 0.0f}, 0.5f };
    coincidence_acceptance("minimal_coincident_opposite_velocity", &C, &D, &c);

    /* 反例：球心完全重合且两球速度相同 -> 相对速度同样退化，必须固定轴 */
    Ball E = { {-2.0f, 1.0f, 0.0f}, {0.7f, -0.7f, 0.0f}, 0.5f };
    Ball F = { {-2.0f, 1.0f, 0.0f}, {0.7f, -0.7f, 0.0f}, 0.5f };
    coincidence_acceptance("counterexample_coincident_equal_velocity", &E, &F, &c);
    check("counterexample_uses_fixed_axis", c.mode == MODE_FIXED_AXIS);
    check("counterexample_normal_axis_x",
          fabsf(c.normal.x - 1.0f) < 1e-6f && c.normal.y == 0.0f && c.normal.z == 0.0f);

    /* 分离的球：无接触，不产生法线 */
    Ball G = { {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, 0.5f };
    Ball H = { {5.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, 0.5f };
    check("separated_no_hit", balls_contact(&G, &H, EPS_T_DEFAULT, &c) == 0);

    /* 临界：d 恰好等于 eps_t，按非退化走几何法线 */
    Ball I = { {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, 0.5f };
    Ball J = { {EPS_T_DEFAULT, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, 0.5f };
    check("boundary_uses_position",
          balls_contact(&I, &J, EPS_T_DEFAULT, &c) == 1 && c.mode == MODE_POSITION);

    /* 非法半径拒绝 */
    Ball K = { {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, 0.0f };
    check("zero_radius_rejected", balls_contact(&K, &A, EPS_T_DEFAULT, &c) == -1);

    /* 深度一致性：重合球深度 = rA+rB */
    balls_contact(&C, &D, EPS_T_DEFAULT, &c);
    check("coincident_depth_full", fabsf(c.depth - 1.0f) < 1e-6f);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
