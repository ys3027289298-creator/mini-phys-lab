
/* naive_repro: 朴素球-球接触生成，无退化处理，用于复现除零。
 * 打开 FE_DIVBYZERO/FE_INVALID 陷阱：0/0 直接 SIGFPE，而不是悄悄变 NaN。 */
#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <fenv.h>

typedef struct { float x, y, z; } V3;
typedef struct { V3 c; V3 v; float r; } Ball;
typedef struct { V3 normal; float depth; int hit; } Contact;

static int contact_naive(const Ball* a, const Ball* b, Contact* o) {
    float dx = b->c.x - a->c.x;
    float dy = b->c.y - a->c.y;
    float dz = b->c.z - a->c.z;
    float d = sqrtf(dx * dx + dy * dy + dz * dz);
    float rsum = a->r + b->r;
    if (d >= rsum) { o->hit = 0; return 0; }
    o->normal.x = dx / d;              /* d==0 时 0/0 */
    o->normal.y = dy / d;
    o->normal.z = dz / d;
    o->depth = rsum - d;
    o->hit = 1;
    return 1;
}

static void run(const char* tag, const Ball* a, const Ball* b) {
    Contact c = {0};
    printf("[%s] running ...\n", tag);
    int h = contact_naive(a, b, &c);
    printf("[%s] hit=%d normal=(%g,%g,%g) depth=%g\n",
           tag, h, c.normal.x, c.normal.y, c.normal.z, c.depth);
}

int main(void) {
    feenableexcept(FE_DIVBYZERO | FE_INVALID);

    /* 原始输入：d=0.02，r=0.5，原点附近，float 尚能分辨 */
    Ball A = { {0.0f, 0.0f, 0.0f}, { 1.0f, 0.0f, 0.0f}, 0.5f };
    Ball B = { {0.02f, 0.0f, 0.0f}, {-1.0f, 0.0f, 0.0f}, 0.5f };
    run("original d=0.02 near origin", &A, &B);

    /* 同一对球，平移到大坐标处：1e6 处 float ULP≈0.0625，0.02 被舍掉，d 直接变 0 */
    Ball A2 = { {1.0e6f, 0.0f, 0.0f}, { 1.0f, 0.0f, 0.0f}, 0.5f };
    Ball B2 = { {1.0e6f + 0.02f, 0.0f, 0.0f}, {-1.0f, 0.0f, 0.0f}, 0.5f };
    run("original d=0.02 at large offset", &A2, &B2);

    /* 最小输入：球心距离为零、半径相等 */
    Ball C = { {3.0f, 4.0f, 5.0f}, { 1.0f, 0.0f, 0.0f}, 0.5f };
    Ball D = { {3.0f, 4.0f, 5.0f}, {-1.0f, 0.0f, 0.0f}, 0.5f };
    run("minimal d=0 equal radii", &C, &D);

    return 0;
}
