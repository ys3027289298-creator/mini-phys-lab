/* integrate.h: the public boundary of the integrate module.

   BodyType is a REQUIRED argument with NO default, and is a struct wrapper
   rather than a bare enum so that a scalar literal cannot stand in for it.
   A call site that does not name the body kind does not compile. */
#ifndef INTEGRATE_H
#define INTEGRATE_H

/* Shared with solve.h via BODY_TYPE_H, so a translation unit may include
   both module headers without a redeclaration. Keep the two copies identical. */
#ifndef BODY_TYPE_H
#define BODY_TYPE_H
typedef enum { BK_DYNAMIC, BK_KINEMATIC } BodyKind;
typedef struct { BodyKind kind; } BodyType;
#define BODY_DYNAMIC   ((BodyType){ BK_DYNAMIC })
#define BODY_KINEMATIC ((BodyType){ BK_KINEMATIC })
#endif

typedef struct { float x, y, px, py, vx, vy; } Body;

/* advance one body by dt.
   Dynamic: gravity integrated, position corrected out of the floor.
   Kinematic: no gravity, no positional correction - its path is authoritative.
   Returns 1 if the body was corrected out of penetration, 0 if not, -1 on bad args. */
int step(BodyType t, Body* b, float dt, float gx, float gy, float floor_y, float radius);

#endif
