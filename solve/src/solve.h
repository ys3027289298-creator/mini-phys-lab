/* solve.h: the public boundary of the solve module.

   BodyType is a REQUIRED argument and has NO default. A call site that does
   not state whether the body is dynamic or kinematic does not compile, so the
   incident class "kinematic body silently solved as dynamic" cannot reach
   the solver.

   The type is a struct wrapper rather than a bare enum on purpose: a scalar
   literal (0, 1) is not convertible to it, so the body kind can only be named
   via BODY_DYNAMIC / BODY_KINEMATIC. A bare enum would accept a literal
   silently - even under gcc/clang with every warning on. */
#ifndef SOLVE_H
#define SOLVE_H

/* Shared with integrate.h via BODY_TYPE_H, so a translation unit may include
   both module headers without a redeclaration. Keep the two copies identical. */
#ifndef BODY_TYPE_H
#define BODY_TYPE_H
typedef enum { BK_DYNAMIC, BK_KINEMATIC } BodyKind;
typedef struct { BodyKind kind; } BodyType;
#define BODY_DYNAMIC   ((BodyType){ BK_DYNAMIC })
#define BODY_KINEMATIC ((BodyType){ BK_KINEMATIC })
#endif

typedef struct { float x, y, z; } V3;
typedef struct { float restitution; float friction; } Material;

/* resolve one contact against a unit normal pointing out of the surface.
   Kinematic bodies have infinite mass: *v is left untouched and 0 is returned. */
int resolve(BodyType t, V3* v, V3 n, const Material* m, float* out_normal_loss);

#endif
