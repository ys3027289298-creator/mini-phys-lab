
/* nan_demo: 不开浮点陷阱时，0/0 产出 NaN，演示零法线如何污染下游 solve。 */
#include <stdio.h>
#include <math.h>

typedef struct { float x, y, z; } V3;

int main(void) {
    /* 最小输入：球心完全重合 */
    V3 da = {0.0f, 0.0f, 0.0f};
    float d = sqrtf(da.x * da.x + da.y * da.y + da.z * da.z);
    V3 n = { da.x / d, da.y / d, da.z / d };
    printf("d=%g  naive normal=(%g,%g,%g)\n", d, n.x, n.y, n.z);
    printf("normal is NaN: %d\n", isnan(n.x));

    /* 下游 solve：法向相对速度 vn = v . n，恢复系数冲量全靠它 */
    V3 v = { -1.0f, 0.0f, 0.0f };
    float vn = v.x * n.x + v.y * n.y + v.z * n.z;
    float new_vn = -vn * 0.5f;
    printf("vn=%g  bounced vn=%g  (vn>0 separating check: %d)\n",
           vn, new_vn, vn > 0.0f);

    /* 单位法线自检也会失效：NaN 比较恒假，防线被静默绕过 */
    float nl = sqrtf(n.x * n.x + n.y * n.y + n.z * n.z);
    printf("|normal|=%g  unit check passes: %d\n", nl, fabsf(nl - 1.0f) < 1e-5f);
    return 0;
}
