/*
 * degenerate: sphere-vs-sphere contact normal when the two centres coincide.
 *
 * The textbook normal n = (cB - cA) / |cB - cA| is a 0/0 form at coincidence.
 * This module defines the degenerate case explicitly with a four-level fallback:
 *
 *   NORMAL_FROM_POSITION  separated/penetrating pair, direction from centres
 *   NORMAL_FROM_VELOCITY  coincident now, limit taken from current rel. velocity
 *   NORMAL_FROM_HISTORY   coincident now and rel. velocity ~ 0, walk back frames
 *   NORMAL_FROM_FIXED_AXIS all else failed; deterministic, but only slides apart
 *
 * History is a per-pair ring buffer of past relative states. Its depth is the
 * only knob the caller has to tune: see experiments/ for the recovery/memory
 * trade-off that picks it.
 */
#ifndef DEGENERATE_H
#define DEGENERATE_H

#include <stddef.h>

#define V3_ZERO ((V3){0.0f, 0.0f, 0.0f})

typedef struct { float x, y, z; } V3;
typedef struct { V3 c; float r; } Sphere;

typedef enum {
    NORMAL_FROM_POSITION   = 0,
    NORMAL_FROM_VELOCITY   = 1,
    NORMAL_FROM_HISTORY    = 2,
    NORMAL_FROM_FIXED_AXIS = 3
} NormalSource;

typedef struct {
    V3     normal;        /* unit vector, points from B toward A (push-out dir) */
    float  depth;         /* rA + rB - |cB - cA|, > 0 means overlap */
    int    hit;           /* 1 while the spheres touch or overlap */
    NormalSource source;  /* which level of the fallback produced n */
    int    history_age;   /* frames walked back; 0 unless source == HISTORY */
} Contact;

/* One past frame of relative motion for a sphere pair. */
typedef struct { V3 rel_pos; V3 rel_vel; } PairFrame;

/* Per-pair ring buffer. capacity is the history depth H. */
typedef struct {
    PairFrame *ring;
    int        capacity;
    int        head;      /* next slot to write */
    int        count;     /* valid frames, saturates at capacity */
} PairHistory;

void  history_init(PairHistory *h, PairFrame *storage, int capacity);
void  history_push(PairHistory *h, V3 rel_pos, V3 rel_vel);

/* Scales with body size so the same code works across radius magnitudes. */
extern float deg_pos_eps_ratio;   /* default 1e-6 */
extern float deg_vel_eps;         /* default 1e-6 m/s */
extern V3    deg_fixed_axis;      /* default +x */

/*
 * Robust contact. Never returns a zero/NaN normal while hit == 1.
 * rel_vel is the current relative velocity (vB - vA) before the impulse.
 */
Contact sphere_contact(const Sphere *a, const Sphere *b,
                       V3 rel_vel, const PairHistory *h);

/*
 * Textbook formula without any degenerate handling. At exact coincidence it
 * performs 0/0: the output normal is NaN (and therefore so is the impulse).
 * Kept to demonstrate the original crash, not for production use.
 */
Contact sphere_contact_naive(const Sphere *a, const Sphere *b);

/* Equal-mass impulse along the contact normal; mirrors solve/ conventions. */
void apply_impulse(V3 *va, V3 *vb, V3 n, float restitution);

/* 1 if any component is NaN or infinite. */
int v3_is_bad(V3 v);

/* Unit length within tolerance, finite components. */
int v3_is_unit(V3 v, float tol);

#endif
