#ifndef PEN_MODEL_H
#define PEN_MODEL_H

/*
 * Shared model behind the penetration spec: projection, camera tiers, and a
 * single 60-second settling run.  Both the acceptance test (src/penetration.c)
 * and the curve dump (tools/sweep.c) include this file, so the numbers the
 * spec is derived from and the numbers the test enforces cannot drift apart.
 */
#include <math.h>
#include "pen_depth.h"

/* ---------------------------------------------------------- projection --- */
/*
 * A penetration of `m` metres at camera distance R covers
 *
 *     px = m * f / R,        f = sensor_px / (2 tan(vfov/2))
 *
 * pixels.  Every threshold in this module is stated in pixels and converted
 * to metres per tier.  f is the only place resolution and FOV enter.
 */
#define PEN_PI 3.14159265358979323846

static inline double pen_focal_px(double sensor_px, double vfov_deg) {
    return sensor_px / (2.0 * tan(vfov_deg * PEN_PI / 360.0));
}

static inline double pen_px_per_m(double R, double f) {
    return f / R;
}

/* --------------------------------------------------------------- tiers --- */
typedef struct {
    double D;      /* horizontal camera distance to the contact, m          */
    double R;      /* slant distance hypot(D, camera height), m             */
    double A_vis;  /* penetration that reaches K_VIS_PX on screen, m        */
    double allow;  /* residual quota chosen for this tier, m                */
} Tier;

typedef struct {
    double radius;
    double plane_y;
    double gravity;     /* magnitude, m/s^2                                */
    double camera_h;    /* camera height above the floor, m                */
    double sensor_px;   /* vertical resolution                             */
    double vfov_deg;    /* vertical field of view                          */
    int    hz;
} World;

static inline void pen_tier(Tier* t, double D, const World* w, double f,
                            double k_vis_px, double k_allow_px) {
    t->D = D;
    t->R = hypot(D, w->camera_h);
    t->A_vis = k_vis_px * t->R / f;
    t->allow = k_allow_px * t->R / f;
}

/* ------------------------------------------------------------- one run --- */
typedef struct {
    long   frames;
    long   visible_frames;     /* rendered penetration > A_vis              */
    long   visible_episodes;   /* contiguous stretches of the above         */
    long   correction_frames;
    double max_physical;       /* worst depth before correction, m          */
    double max_render;         /* worst rendered depth, m                   */
    double max_render_px;      /* ... in pixels at this tier's R            */
    double final_render;       /* rendered depth on the last frame, m       */
} Run;

/*
 * Per-frame order, matching ../integrate (correction happens inside step,
 * before the frame is presented):
 *
 *   v += g*dt;  y += v*dt;                        integrate
 *   h  = pen_depth(y)                              physical penetration
 *   if (h > allow) { y += beta*h; }                position correction
 *   hr = pen_depth(y)                              what the camera sees
 *   visible if hr*A_vis^{-1} > 1                   screen test
 *   if (h > 0 && v < 0) v = -e*v                   bounce
 *
 * beta = 1 snaps the ball back to the surface on the correction frame, so the
 * residual a viewer can see is exactly the un-corrected growth between
 * corrections.  `allow` is the deadband: penetration at or below it is left
 * alone, which is what lets a run skip correction frames.  With beta = 1 the
 * rendered depth never exceeds `allow`, so the only question the test has to
 * answer is whether `allow` maps to more than one screen pixel at distance R.
 */
static inline Run pen_run(const Tier* t, double allow, double beta, double e,
                          double drop_m, const World* w) {
    Run r;
    r.frames = r.visible_frames = r.visible_episodes = r.correction_frames = 0;
    r.max_physical = r.max_render = r.max_render_px = r.final_render = 0.0;

    double dt = 1.0 / w->hz;
    long n = (long)(w->hz * 60.0);          /* 60 s of settling per tier   */
    double y = w->plane_y + w->radius + drop_m;
    double v = 0.0;
    int was_visible = 0;
    r.frames = n;

    for (long i = 0; i < n; i++) {
        v += -w->gravity * dt;
        y += v * dt;

        double h = pen_depth(y, w->radius, w->plane_y);
        if (h > r.max_physical) r.max_physical = h;

        /* Snap correction; the clamp keeps beta>1 an over-relaxation rather
         * than a teleport through the surface. */
        if (h > allow) {
            double c = beta * h;
            if (c > h) c = h;
            y += c;
            r.correction_frames++;
        }

        double hr = pen_depth(y, w->radius, w->plane_y);
        if (hr > r.max_render) r.max_render = hr;
        double px = hr * pen_px_per_m(t->R, pen_focal_px(w->sensor_px, w->vfov_deg));
        if (px > r.max_render_px) r.max_render_px = px;
        r.final_render = hr;

        int vis = (hr > t->A_vis);
        if (vis) r.visible_frames++;
        if (vis && !was_visible) r.visible_episodes++;
        was_visible = vis;

        if (h > 0.0 && v < 0.0) v = -v * e;
    }
    return r;
}

#endif
