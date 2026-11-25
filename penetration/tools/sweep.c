/*
 * sweep: dump the two curves the quota was chosen from.
 *
 * For a grid of residual quotas (expressed in screen pixels at each tier,
 * because that is the scale that matters), report
 *   - visible episode count   (how often a viewer can see the ball sunk)
 *   - correction frame count  (what the residual policy costs)
 *
 * The chosen quota is the loosest point whose visible count is still zero,
 * which on the 3 m tier is pinned to the < 1 px hard limit regardless of cost.
 *
 * Output: CSV on stdout, comments on stderr.
 */
#include <stdio.h>
#include "pen_model.h"

#define CAMERA_H    1.50
#define SENSOR_PX   1080.0
#define VFOV_DEG    60.0
#define RADIUS      0.50
#define PLANE_Y     0.0
#define GRAVITY     9.8
#define DROP        6.00
#define RESTITUTION 0.20
#define HZ          60

static const double TIER_D[3] = { 3.0, 10.0, 30.0 };
static const char*  TIER_NAME[3] = { "3m", "10m", "30m" };

int main(void) {
    World w = { RADIUS, PLANE_Y, GRAVITY, CAMERA_H, SENSOR_PX, VFOV_DEG, HZ };
    printf("tier,D_m,R_m,quota_px,quota_mm,quota_frac_radius,"
           "visible_episodes,visible_frames,correction_frames,peak_px\n");

    for (int i = 0; i < 3; i++) {
        double R = hypot(TIER_D[i], w.camera_h);
        double f = pen_focal_px(w.sensor_px, w.vfov_deg);
        /* quota from 0.05 px to 4.0 px of on-screen penetration */
        for (int k = 1; k <= 80; k++) {
            double qpx = 0.05 * k;
            double allow = qpx * R / f;             /* px -> metres */
            if (allow > RADIUS) continue;
            Tier t;
            pen_tier(&t, TIER_D[i], &w, f, 1.0, qpx);
            Run r = pen_run(&t, allow, 1.00, RESTITUTION, DROP, &w);
            printf("%s,%.1f,%.4f,%.2f,%.3f,%.4f,%ld,%ld,%ld,%.3f\n",
                   TIER_NAME[i], TIER_D[i], R, qpx, allow * 1e3,
                   allow / RADIUS, r.visible_episodes, r.visible_frames,
                   r.correction_frames, r.max_render_px);
        }
    }
    return 0;
}
