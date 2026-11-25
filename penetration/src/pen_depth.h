#ifndef PEN_DEPTH_H
#define PEN_DEPTH_H

/*
 * Penetration depth of a sphere against the half-space y >= plane_y.
 *
 *   h = plane_y + radius - centre_y
 *
 * h > 0  -> the sphere is below the surface by h metres (penetrating).
 * h == 0 -> resting contact.
 * h < 0  -> clear of the surface.
 *
 * This is the only definition used by the module; the checks and the
 * numeric spec both go through it, so they cannot drift apart.
 */
static inline double pen_depth(double centre_y, double radius, double plane_y) {
    return (plane_y + radius) - centre_y;
}

#endif
