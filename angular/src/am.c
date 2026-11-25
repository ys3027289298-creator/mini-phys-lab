
/* am: rigid-body normal impulse with explicit angular-momentum bookkeeping.
 * See include/am.h for the sign convention and the conservation inequality. */
#include "am.h"

#include <string.h>

am_vec3 am_v(double x, double y, double z) { am_vec3 r = { x, y, z }; return r; }
am_vec3 am_add(am_vec3 a, am_vec3 b) { return am_v(a.x + b.x, a.y + b.y, a.z + b.z); }
am_vec3 am_sub(am_vec3 a, am_vec3 b) { return am_v(a.x - b.x, a.y - b.y, a.z - b.z); }
am_vec3 am_scale(am_vec3 a, double s) { return am_v(a.x * s, a.y * s, a.z * s); }
double  am_dot(am_vec3 a, am_vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
double  am_len2(am_vec3 a) { return am_dot(a, a); }
double  am_len(am_vec3 a) { return sqrt(am_len2(a)); }

/* right-handed cross product: ex x ey = ez */
am_vec3 am_cross(am_vec3 a, am_vec3 b) {
    return am_v(a.y * b.z - a.z * b.y,
                a.z * b.x - a.x * b.z,
                a.x * b.y - a.y * b.x);
}

static double inv_mass(const am_body* b) { return (b && b->m > 0.0) ? 1.0 / b->m : 0.0; }
static double inv_inertia(const am_body* b) {
    return (b && b->m > 0.0 && b->I > 0.0) ? 1.0 / b->I : 0.0;
}

static am_vec3 contact_vel(const am_body* b, am_vec3 r) {
    /* u = v + omega x r */
    return am_add(b->v, am_cross(b->w, r));
}

am_vec3 am_total_L(const am_body* A, const am_body* B, am_vec3 o) {
    am_vec3 L = am_v(0, 0, 0);
    const am_body* bs[2] = { A, B };
    int i;
    for (i = 0; i < 2; ++i) {
        const am_body* b = bs[i];
        if (!b) continue;
        am_vec3 d = am_sub(b->x, o);
        /* orbital: d x (m v); static bodies carry no finite momentum */
        if (b->m > 0.0) L = am_add(L, am_cross(d, am_scale(b->v, b->m)));
        /* spin: I omega */
        if (b->I > 0.0) L = am_add(L, am_scale(b->w, b->I));
    }
    return L;
}

double am_gross_scale(const am_body* A, const am_body* B, am_vec3 o) {
    double s = 0.0;
    const am_body* bs[2] = { A, B };
    int i;
    for (i = 0; i < 2; ++i) {
        const am_body* b = bs[i];
        if (!b) continue;
        am_vec3 d = am_sub(b->x, o);
        if (b->m > 0.0) s += am_len(am_cross(d, am_scale(b->v, b->m)));
        if (b->I > 0.0) s += fabs(b->I * am_len(b->w));
    }
    return s;
}

am_snapshot am_snapshot_make(const am_body* A, const am_body* B) {
    am_snapshot s;
    double mA = (A && A->m > 0.0) ? A->m : 0.0;
    double mB = (B && B->m > 0.0) ? B->m : 0.0;
    double mt = mA + mB;
    if (mt > 0.0) {
        am_vec3 num = am_add(am_scale(A->x, mA), am_scale(B->x, mB));
        s.origin = am_scale(num, 1.0 / mt);
    } else {
        s.origin = am_v(0, 0, 0);
    }
    s.L_before = am_total_L(A, B, s.origin);
    s.total_before = am_len(s.L_before);
    s.gross = am_gross_scale(A, B, s.origin);
    return s;
}

/* ---- float32 tier: same math, single precision, one cross per body ------ */
static am_solution solve_apply_f32(am_body* A, am_body* B, const am_contact* c) {
    am_solution sol;
    memset(&sol, 0, sizeof(sol));
    sol.status = 0;

    float nX = (float)c->n.x, nY = (float)c->n.y, nZ = (float)c->n.z;
    float nl = sqrtf(nX * nX + nY * nY + nZ * nZ);
    if (nl < 1e-20f) { sol.status = -2; return sol; }
    nX /= nl; nY /= nl; nZ /= nl;

    float rAx = (float)(c->p.x - A->x.x), rAy = (float)(c->p.y - A->x.y), rAz = (float)(c->p.z - A->x.z);
    float rBx = (float)(c->p.x - B->x.x), rBy = (float)(c->p.y - B->x.y), rBz = (float)(c->p.z - B->x.z);

    /* cA = rA x n (z components written fully for clarity) */
    float cAx = rAy * nZ - rAz * nY, cAy = rAz * nX - rAx * nZ, cAz = rAx * nY - rAy * nX;
    float cBx = rBy * nZ - rBz * nY, cBy = rBz * nX - rBx * nZ, cBz = rBx * nY - rBy * nX;
    float eccA = sqrtf(cAx * cAx + cAy * cAy + cAz * cAz);
    float eccB = sqrtf(cBx * cBx + cBy * cBy + cBz * cBz);
    sol.eccA = eccA; sol.eccB = eccB;
    sol.rA = am_v(rAx, rAy, rAz); sol.rB = am_v(rBx, rBy, rBz);

    /* contact velocities (float) */
    float wAx = (float)A->w.x, wAy = (float)A->w.y, wAz = (float)A->w.z;
    float uAx = (float)A->v.x + (wAy * rAz - wAz * rAy);
    float uAy = (float)A->v.y + (wAz * rAx - wAx * rAz);
    float uAz = (float)A->v.z + (wAx * rAy - wAy * rAx);
    float wBx = (float)B->w.x, wBy = (float)B->w.y, wBz = (float)B->w.z;
    float uBx = (float)B->v.x + (wBy * rBz - wBz * rBy);
    float uBy = (float)B->v.y + (wBz * rBx - wBx * rBz);
    float uBz = (float)B->v.z + (wBx * rBy - wBy * rBx);

    float dux = uBx - uAx, duy = uBy - uAy, duz = uBz - uAz;
    float vn = dux * nX + duy * nY + duz * nZ;
    sol.vn = vn;
    if (vn >= 0.0f) { sol.status = 0; return sol; }

    float imA = (A->m > 0.0) ? 1.0f / (float)A->m : 0.0f;
    float imB = (B->m > 0.0) ? 1.0f / (float)B->m : 0.0f;
    float iiA = (A->m > 0.0 && A->I > 0.0) ? 1.0f / (float)A->I : 0.0f;
    float iiB = (B->m > 0.0 && B->I > 0.0) ? 1.0f / (float)B->I : 0.0f;
    float k = imA + imB
            + (cAx * cAx + cAy * cAy + cAz * cAz) * iiA
            + (cBx * cBx + cBy * cBy + cBz * cBz) * iiB;
    if (!(k > 0.0f)) { sol.status = -3; return sol; }
    float e = (float)c->e;
    float j = -(1.0f + e) * vn / k;
    sol.j = j; sol.k_eff = k; sol.status = 1;

    /* impulse on A = -j n, on B = +j n */
    if (imA > 0.0f) {
        A->v.x += (double)(-j * nX * imA);
        A->v.y += (double)(-j * nY * imA);
        A->v.z += (double)(-j * nZ * imA);
    }
    if (iiA > 0.0f) { /* I dw = r x J_A, one explicit lever-arm cross */
        float jAx = -j * nX, jAy = -j * nY, jAz = -j * nZ;
        float tx = rAy * jAz - rAz * jAy;
        float ty = rAz * jAx - rAx * jAz;
        float tz = rAx * jAy - rAy * jAx;
        A->w.x += (double)(tx * iiA);
        A->w.y += (double)(ty * iiA);
        A->w.z += (double)(tz * iiA);
    }
    if (imB > 0.0f) {
        B->v.x += (double)(j * nX * imB);
        B->v.y += (double)(j * nY * imB);
        B->v.z += (double)(j * nZ * imB);
    }
    if (iiB > 0.0f) { /* I dw = r x J_B */
        float jBx = j * nX, jBy = j * nY, jBz = j * nZ;
        float tx = rBy * jBz - rBz * jBy;
        float ty = rBz * jBx - rBx * jBz;
        float tz = rBx * jBy - rBy * jBx;
        B->w.x += (double)(tx * iiB);
        B->w.y += (double)(ty * iiB);
        B->w.z += (double)(tz * iiB);
    }
    return sol;
}

am_solution am_solve_apply(am_body* A, am_body* B, const am_contact* c,
                           am_tier tier, am_vec3 errA, am_vec3 errB) {
    am_solution sol;
    memset(&sol, 0, sizeof(sol));
    sol.status = 0;
    if (!A || !B || !c) { sol.status = -1; return sol; }
    if (c->e < 0.0 || c->e > 1.0) { sol.status = -1; return sol; }
    if (tier == AM_TIER_F32) return solve_apply_f32(A, B, c);

    /* ---- double path (also the coarse particle path) ---- */
    am_vec3 n = c->n;
    double nl = am_len(n);
    if (nl < 1e-18) { sol.status = -2; return sol; }
    n = am_scale(n, 1.0 / nl);

    am_vec3 rA = am_sub(c->p, A->x);
    am_vec3 rB = am_sub(c->p, B->x);
    sol.rA = rA; sol.rB = rB;

    int rotate = (tier != AM_TIER_COARSE);
    /* effective-mass arm uses the TRUE geometry; injected error belongs to
     * the angular-impulse arm (the precision being audited). */
    am_vec3 cA = am_cross(rA, n);
    am_vec3 cB = am_cross(rB, n);
    sol.eccA = am_len(cA);
    sol.eccB = am_len(cB);

    /* The coarse particle tier ignores rotation entirely: contact velocity is
     * just v and no lever-arm cross product is evaluated. */
    am_vec3 uA = rotate ? contact_vel(A, rA) : A->v;
    am_vec3 uB = rotate ? contact_vel(B, rB) : B->v;
    double vn = am_dot(am_sub(uB, uA), n);
    sol.vn = vn;
    if (vn >= 0.0) { sol.status = 0; return sol; } /* separating */

    double imA = inv_mass(A), imB = inv_mass(B);
    double iiA = rotate ? inv_inertia(A) : 0.0;
    double iiB = rotate ? inv_inertia(B) : 0.0;

    double k = imA + imB;
    if (rotate) {
        k += am_len2(cA) * iiA + am_len2(cB) * iiB;
    }
    if (!(k > 0.0)) { sol.status = -3; return sol; }

    double j = -(1.0 + c->e) * vn / k;
    am_vec3 JA = am_scale(n, -j);   /* impulse ON A */
    am_vec3 JB = am_scale(n, +j);   /* impulse ON B */
    sol.j = j; sol.k_eff = k; sol.status = 1;

    if (imA > 0.0) A->v = am_add(A->v, am_scale(JA, imA));
    if (imB > 0.0) B->v = am_add(B->v, am_scale(JB, imB));

    if (rotate) {
        if (iiA > 0.0) {
            am_vec3 armA = am_add(rA, errA);
            A->w = am_add(A->w, am_scale(am_cross(armA, JA), iiA));
        }
        if (iiB > 0.0) {
            am_vec3 armB = am_add(rB, errB);
            B->w = am_add(B->w, am_scale(am_cross(armB, JB), iiB));
        }
    }
    return sol;
}

am_report am_check(const am_snapshot* snap,
                   const am_body* A1, const am_body* B1,
                   double eccA, double eccB, double tol) {
    am_report r;
    memset(&r, 0, sizeof(r));
    r.tol = (tol > 0.0) ? tol : AM_TOL_ONE_PER_MILLE;
    r.L_before = snap->L_before;
    r.total_before = snap->total_before;
    r.gross = snap->gross;

    am_vec3 L_after = am_total_L(A1, B1, snap->origin);
    r.L_after = L_after;
    r.dL = am_sub(L_after, snap->L_before);
    r.dL_abs = am_len(r.dL);

    r.rel_total = r.dL_abs / (r.total_before + 1e-300);
    r.rel_gross = r.dL_abs / (r.gross + 1e-300);
    r.near_cancel = (r.total_before < 0.1 * r.gross);

    r.eccA = eccA; r.eccB = eccB;
    r.ecc_max = (eccA > eccB) ? eccA : eccB;
    r.pass = (r.rel_gross <= r.tol);
    return r;
}

void am_canonical(am_body* A, am_body* B, am_contact* c) {
    /* Hand-checkable eccentric pair. Contact at origin, n = +x (A left -> B).
     * rA = (+0.5,+0.30,0), rB = (-0.5,-0.15,0); eccA=0.30, eccB=0.15.
     * mA=1, I_A=0.1 ; mB=2, I_B=0.2 ; vn = -10 ; e=0.5.
     * COM = (1/6,0,0), L_before = (0,0,3.0). */
    memset(A, 0, sizeof(*A));
    memset(B, 0, sizeof(*B));
    A->m = 1.0; A->I = 0.1;
    A->x = am_v(-0.5, -0.30, 0.0);
    A->v = am_v(8.0, 0.0, 0.0);
    A->w = am_v(0, 0, 0);

    B->m = 2.0; B->I = 0.2;
    B->x = am_v(0.5, 0.15, 0.0);
    B->v = am_v(-2.0, 0.0, 0.0);
    B->w = am_v(0, 0, 0);

    c->p = am_v(0, 0, 0);
    c->n = am_v(1, 0, 0);
    c->e = 0.5;
}
