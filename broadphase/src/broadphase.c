
/* broadphase: uniform-grid broad phase for sphere-sphere collision in 2D.
 *
 * Broad phase keeps every pair whose grid cells touch (the spatial-hash form
 * of AABB culling); the narrow phase then runs the exact distance test
 * d <= r_i + r_j.  AABB overlap is only necessary, never sufficient: two
 * boxes can overlap at a corner while the spheres are farther than r_i+r_j,
 * so every surviving pair still needs a narrow-phase check.
 *
 * Modes:
 *   ./broadphase test            self-check (deterministic, no arguments)
 *   ./broadphase sweep <file>    scan cell size h for three radius tiers,
 *                                write one CSV row per (tier, h)
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef struct { double x, y, r; } Body;

/* ---------- deterministic RNG (xorshift32) ---------- */
static unsigned long rng_state;
static double rng_next(void) {
    unsigned long x = rng_state;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    rng_state = x;
    return (double)(x & 0xFFFFFFu) / (double)0x1000000u; /* [0,1) */
}

static void make_scene(Body* b, int n, double L, double r, unsigned long seed) {
    rng_state = seed;
    for (int i = 0; i < n; i++) {
        b[i].x = rng_next() * L;
        b[i].y = rng_next() * L;
        b[i].r = r;
    }
}

/* ---------- exhaustive reference: every pair, exact distance test ---------- */
static long brute_true_pairs(const Body* b, int n) {
    long true_pairs = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double dx = b[i].x - b[j].x, dy = b[i].y - b[j].y;
            double rr = b[i].r + b[j].r;
            if (dx * dx + dy * dy <= rr * rr) true_pairs++;
        }
    }
    return true_pairs;
}

typedef struct {
    long candidates;   /* pairs surviving broad phase -> narrow-phase calls */
    long true_pairs;   /* candidates that actually collide */
    long misses;       /* real collisions rejected by broad phase */
    long false_pos;    /* candidates that survive but do not collide */
} BpResult;

/*
 * Uniform grid.  Each body is keyed by the cell its centre falls in; a body
 * only checks the 3x3 neighbourhood of its own cell.  Pair (i,j), i<j, is
 * counted once, when i is visited and j is one of its neighbours.
 *
 * Zero-miss condition: colliding centres are at most r_i+r_j apart, which is
 * < 2*h whenever h >= max(r_i+r_j); their cell indices then differ by at
 * most 1 in each axis.  With h < 2*r_max some colliding pairs land in cells
 * two apart and are dropped -> misses.
 */
static BpResult grid_broadphase(const Body* b, int n, double L, double h) {
    BpResult res = {0, 0, 0, 0};
    if (h <= 0.0) { res.misses = brute_true_pairs(b, n); return res; }

    long g = (long)floor(L / h);    if (g < 1) g = 1;
    double cell = L / (double)g;           /* actual cell size, divides L */

    long* head = (long*)malloc((size_t)g * g * sizeof(long));
    long* next = (long*)malloc((size_t)n * sizeof(long));
    if (!head || !next) { free(head); free(next); exit(3); }
    for (long c = 0; c < g * g; c++) head[c] = -1;

    for (int i = 0; i < n; i++) {
        long cx = (long)(b[i].x / cell); if (cx >= g) cx = g - 1;
        long cy = (long)(b[i].y / cell); if (cy >= g) cy = g - 1;
        long c = cy * g + cx;
        next[i] = head[c];
        head[c] = i;
    }

    long total_true = brute_true_pairs(b, n);

    for (int i = 0; i < n; i++) {
        long cx = (long)(b[i].x / cell); if (cx >= g) cx = g - 1;
        long cy = (long)(b[i].y / cell); if (cy >= g) cy = g - 1;
        for (long oy = -1; oy <= 1; oy++) {
            long ny = cy + oy;
            if (ny < 0 || ny >= g) continue;
            for (long ox = -1; ox <= 1; ox++) {
                long nx = cx + ox;
                if (nx < 0 || nx >= g) continue;
                for (long j = head[ny * g + nx]; j != -1; j = next[j]) {
                    if (j <= i) continue;                 /* each pair once */
                    res.candidates++;
                    double dx = b[i].x - b[j].x, dy = b[i].y - b[j].y;
                    double rr = b[i].r + b[j].r;
                    if (dx * dx + dy * dy <= rr * rr) res.true_pairs++;
                }
            }
        }
    }

    res.misses = total_true - res.true_pairs;
    res.false_pos = res.candidates - res.true_pairs;
    free(head);
    free(next);
    return res;
}

/* ---------- self-check ---------- */
static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(int argc, char** argv) {
    enum { N = 400 };
    const double L = 100.0;
    const double radii[3] = {0.5, 1.0, 2.0};
    const unsigned long seeds[3] = {0x12345678u, 0x9abcdef0u, 0x55aa5aa5u};

    if (argc >= 2 && strcmp(argv[1], "test") == 0) {
        /* L=10 so h divides L exactly and cell indices are hand-checkable. */
        const double TL = 10.0;
        /* centres 0.2 apart (< 2r): a colliding pair.  At h=2r=2 the centres
         * land in adjacent cells (2 and 3) -> broad phase keeps it. */
        Body cross[2] = { {4.9, 5.0, 1.0}, {5.1, 5.0, 1.0} };
        BpResult r1 = grid_broadphase(cross, 2, TL, 2.0);
        check("boundary_pair_kept_at_h_eq_2r", r1.candidates == 1 && r1.true_pairs == 1 && r1.misses == 0);

        /* colliding pair (centres 1.2 apart < 2r) whose cells are two apart
         * at h=r=1: cell 4 and cell 6, outside the 3x3 window -> dropped. */
        Body gap[2] = { {4.9, 5.0, 1.0}, {6.1, 5.0, 1.0} };
        BpResult r2 = grid_broadphase(gap, 2, TL, 1.0);
        check("boundary_pair_missed_below_h_eq_2r", r2.candidates == 0 && r2.misses == 1);

        /* far pair is never a candidate regardless of h */
        Body far[2] = { {1.0, 1.0, 1.0}, {8.0, 8.0, 1.0} };
        BpResult r3 = grid_broadphase(far, 2, TL, 2.0);
        check("far_pair_culled", r3.candidates == 0 && r3.true_pairs == 0);

        /* diagonally adjacent 1-unit cells: cell AABBs meet at a corner but
         * the sphere centres are ~0.28 apart, farther than r_i+r_j=0.2 ->
         * candidate that the narrow-phase distance test must reject */
        Body corner[2] = { {4.9, 4.9, 0.1}, {5.1, 5.1, 0.1} };
        BpResult r4 = grid_broadphase(corner, 2, TL, 1.0);
        check("corner_overlap_is_false_positive",
              r4.candidates == 1 && r4.true_pairs == 0 && r4.false_pos == 1);

        /* full random scenes: brute force agrees with broad phase at safe h */
        Body* scene = (Body*)malloc((size_t)N * sizeof(Body));
        for (int t = 0; t < 3; t++) {
            make_scene(scene, N, L, radii[t], seeds[t]);
            long truth = brute_true_pairs(scene, N);
            BpResult safe = grid_broadphase(scene, N, L, 2.0 * radii[t]);
            char nm[96];
            snprintf(nm, sizeof(nm), "tier_r%.1f_zero_miss", radii[t]);
            check(nm, safe.misses == 0 && safe.true_pairs == truth && truth > 0);
        }
        free(scene);

        printf("all_ok=%d\n", failures == 0);
        return failures == 0 ? 0 : 1;
    }

    if (argc >= 3 && strcmp(argv[1], "sweep") == 0) {
        FILE* f = fopen(argv[2], "w");
        if (!f) { perror(argv[2]); return 2; }
        fprintf(f, "radius,h,cells_per_axis,candidates,true_pairs,misses,false_pos,"
                   "total_pairs,narrow_frac,reduction_vs_brute,theory_frac\n");

        Body* scene = (Body*)malloc((size_t)N * sizeof(Body));
        long total_pairs = (long)N * (N - 1) / 2;

        for (int t = 0; t < 3; t++) {
            make_scene(scene, N, L, radii[t], seeds[t]);
            long truth = brute_true_pairs(scene, N);
            (void)truth;
            /* geometric h sweep: sqrt(2) steps, spanning both sides of h=2r */
            for (int k = 0; k <= 18; k++) {
                double h = 0.25 * pow(sqrt(2.0), (double)k);
                BpResult r = grid_broadphase(scene, N, L, h);
                long g = (long)fmax(1.0, floor(L / h));
                double frac = (double)r.candidates / (double)total_pairs;
                double theory = 9.0 * (h / L) * (h / L);
                if (theory > 1.0) theory = 1.0;
                fprintf(f, "%.3f,%.6f,%ld,%ld,%ld,%ld,%ld,%ld,%.8f,%.8f,%.8f\n",
                        radii[t], h, g, r.candidates, r.true_pairs, r.misses,
                        r.false_pos, total_pairs, frac, 1.0 - frac, theory);
            }
        }
        fclose(f);
        free(scene);
        printf("sweep written to %s (N=%d L=%.0f total_pairs=%ld)\n",
               argv[2], N, L, total_pairs);
        return 0;
    }

    fprintf(stderr, "usage: %s test | sweep <csv-file>\n", argv[0]);
    return 2;
}
