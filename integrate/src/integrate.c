
/* integrate: semi-implicit euler with position correction against penetration. */
#include <stdio.h>
#include <string.h>
#include <math.h>

typedef struct { float x, y, px, py, vx, vy; } Body;

static int step(Body* b, float dt, float gx, float gy, float floor_y, float radius) {
    if (dt < 0.0f || radius <= 0.0f) return -1;
    b->px = b->x; b->py = b->y;
    b->vx += gx * dt;
    b->vy += gy * dt;
    b->x += b->vx * dt;
    b->y += b->vy * dt;
    /* floor at y = floor_y, body bottom is y - radius */
    float pen = (floor_y + radius) - b->y;
    if (pen > 0.0f) {
        b->y += pen;              /* positional correction */
        if (b->vy < 0.0f) b->vy = 0.0f;
    }
    return pen > 0.0f ? 1 : 0;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Body b; memset(&b, 0, sizeof(b));
    b.y = 10.0f;
    int hit = step(&b, 0.1f, 0, -9.8f, 0.0f, 0.5f);
    check("falling_no_hit", hit == 0 && b.y < 10.0f && b.y > 0.0f);

    /* drive it into the floor */
    for (int i = 0; i < 100; i++) step(&b, 0.1f, 0, -9.8f, 0.0f, 0.5f);
    check("lands_on_floor", fabsf(b.y - 0.5f) < 1e-3f);
    check("vertical_stopped", fabsf(b.vy) < 1e-6f);

    /* never sinks below floor across a long run */
    Body c; memset(&c, 0, sizeof(c));
    c.y = 3.0f; c.vx = 5.0f;
    float min_y = 1e9f;
    for (int i = 0; i < 500; i++) {
        step(&c, 0.05f, 0, -9.8f, 0.0f, 0.5f);
        if (c.y < min_y) min_y = c.y;
    }
    check("never_sinks", min_y >= 0.5f - 1e-4f);

    Body d; memset(&d, 0, sizeof(d));
    check("bad_dt", step(&d, -1.0f, 0, 0, 0, 0.5f) == -1);
    check("bad_radius", step(&d, 0.1f, 0, 0, 0, 0.0f) == -1);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
