/* Fault-injection adapter for integrate. Includes the unmodified module so the
 * static step() is reachable. One case per process for clean gcov. */
#define main integrate_module_main
#include "../../integrate/src/integrate.c"
#undef main
#include <stdlib.h>

static volatile int g_sink;

static void jf(float v) {
    if (isnan(v)) fputs("\"NaN\"", stdout);
    else if (isinf(v)) fputs(v > 0 ? "\"Inf\"" : "\"-Inf\"", stdout);
    else printf("%.9g", v);
}
static void jbody(const Body* b) {
    printf("{\"x\":"); jf(b->x); printf(",\"y\":"); jf(b->y);
    printf(",\"px\":"); jf(b->px); printf(",\"py\":"); jf(b->py);
    printf(",\"vx\":"); jf(b->vx); printf(",\"vy\":"); jf(b->vy); printf("}");
}

/* ---- NUM: infinite dt penetrates the dt-guard ------------------------- */
static void case_N3(int emit) {
    Body b; memset(&b, 0, sizeof b); b.y = 10.0f;
    Body before = b;
    int ret = step(&b, INFINITY, 0.0f, -9.8f, 0.0f, 0.5f);
    g_sink += ret;
    if (!emit) return;
    printf("{\"module\":\"integrate\",\"case\":\"N3\",\"ret\":%d,\"before\":", ret);
    jbody(&before); printf(",\"after\":"); jbody(&b); printf("}\n");
}

/* ---- STATE: a rejected step must not corrupt / leave mixed prev-state ---
 * Fair construction: reach a clean, self-consistent state with one valid step
 * (px,py then legally describe "before the last successful step"), snapshot it,
 * issue a rejected step (radius 0), and compare. A rejected step must be a
 * no-op: x,y,vx,vy and px,py all unchanged. */
static void case_S3(int emit) {
    Body b; memset(&b, 0, sizeof b);
    b.y = 5.0f;
    step(&b, 0.1f, 0.0f, -9.8f, 0.0f, 0.5f);   /* valid step: establishes prev-state */
    Body before = b;
    int ret = step(&b, 0.1f, 0.0f, -9.8f, 0.0f, 0.0f);  /* radius 0 -> -1 */
    g_sink += ret;
    if (!emit) return;
    printf("{\"module\":\"integrate\",\"case\":\"S3\",\"ret\":%d,\"before\":", ret);
    jbody(&before); printf(",\"after\":"); jbody(&b); printf("}\n");
}

/* ---- TIMING: permute independent bodies -------------------------------- */
typedef struct { float y; float gx, gy; } IB;
static const IB t3_bodies[3] = {
    { 10.0f, 0.0f, -9.8f },   /* fall from high      */
    {  0.2f, 0.0f, -9.8f },   /* start below floor   */
    {  3.0f, 5.0f,  0.0f  },  /* sideways drift      */
};
static void t3_run(const int* order, int out_ret[3][2], Body out_b[3]) {
    Body st[3];
    for (int i = 0; i < 3; i++) { memset(&st[i], 0, sizeof st[i]); st[i].y = t3_bodies[i].y; }
    for (int k = 0; k < 3; k++) {
        int id = order[k];
        for (int s = 0; s < 2; s++) {
            out_ret[id][s] = step(&st[id], 0.1f, t3_bodies[id].gx, t3_bodies[id].gy, 0.0f, 0.5f);
        }
    }
    for (int i = 0; i < 3; i++) out_b[i] = st[i];
}
static void case_T3(int emit) {
    static const int ord[3]  = {0,1,2};
    static const int perm[3] = {2,0,1};
    int oret[3][2], pret[3][2]; Body ob[3], pb[3];
    t3_run(ord, oret, ob);
    t3_run(perm, pret, pb);
    float maxdiff = 0.0f;
    for (int i = 0; i < 3; i++) {
        float d = fabsf(ob[i].x - pb[i].x);
        d = fmaxf(d, fabsf(ob[i].y - pb[i].y));
        d = fmaxf(d, fabsf(ob[i].px - pb[i].px));
        d = fmaxf(d, fabsf(ob[i].py - pb[i].py));
        d = fmaxf(d, fabsf(ob[i].vx - pb[i].vx));
        d = fmaxf(d, fabsf(ob[i].vy - pb[i].vy));
        if (oret[i][0] != pret[i][0] || oret[i][1] != pret[i][1]) d = INFINITY;
        maxdiff = fmaxf(maxdiff, d);
        g_sink += oret[i][0] + oret[i][1] + pret[i][0] + pret[i][1];
    }
    if (!emit) return;
    printf("{\"module\":\"integrate\",\"case\":\"T3\",\"maxdiff\":"); jf(maxdiff);
    printf(",\"ordered\":[");
    for (int i = 0; i < 3; i++) { if (i) putchar(','); printf("{\"id\":%d,\"body\":", i); jbody(&ob[i]); putchar('}'); }
    printf("],\"permuted\":[");
    for (int i = 0; i < 3; i++) { if (i) putchar(','); printf("{\"id\":%d,\"body\":", i); jbody(&pb[i]); putchar('}'); }
    printf("]}\n");
}

/* ---- extra candidate injection points (coverage experiment) ----------- */
typedef struct { const char* id; float y, vx; float dt, gx, gy, floor_y, r; int n; } IProbe;
static const IProbe iprobes[] = {
    {"R_fall",    10.0f, 0, 0.1f, 0,-9.8f, 0, 0.5f, 1},
    {"R_land",     0.4f, 0, 0.1f, 0,-9.8f, 0, 0.5f, 3},
    {"R_pen",      0.2f, 0, 0.1f, 0, 0.0f, 0, 0.5f, 1},
    {"R_side",     3.0f, 5, 0.05f,0, 0.0f, 0, 0.5f, 20},
    {"R_baddt",    5.0f, 0, -1.0f,0, 0.0f, 0, 0.5f, 1},
    {"R_badr",     5.0f, 0, 0.1f, 0, 0.0f, 0, 0.0f, 1},
    {"R_nandt",   10.0f, 0, 0.0f/0.0f,0,-9.8f,0,0.5f, 1},
    {"R_infdt",   10.0f, 0, INFINITY,0,-9.8f,0,0.5f, 1},
    {"R_nangx",    1.0f, 0, 0.1f, 0.0f/0.0f,0,0,0.5f, 1},
    {"R_upward",   0.6f, 0, 0.1f, 0, 5.0f, 0, 0.5f, 1},
    {"R_longfall",50.0f, 0, 0.3f, 0,-9.8f, 0, 0.5f, 60},
};
#define NR ((int)(sizeof(iprobes)/sizeof(iprobes[0])))

static void case_rprobe(int idx, int emit) {
    const IProbe* p = &iprobes[idx];
    Body b; memset(&b, 0, sizeof b); b.y = p->y; b.vx = p->vx;
    int ret = 0;
    for (int i = 0; i < p->n; i++) ret = step(&b, p->dt, p->gx, p->gy, p->floor_y, p->r);
    g_sink += ret;
    if (!emit) return;
    printf("{\"module\":\"integrate\",\"case\":\"%s\",\"ret\":%d,\"after\":", p->id, ret);
    jbody(&b); printf("}\n");
}

static int find_r(const char* id) {
    for (int i = 0; i < NR; i++) if (!strcmp(iprobes[i].id, id)) return i;
    return -1;
}
static void run_case(const char* id, int emit) {
    if (!strcmp(id, "N3")) return case_N3(emit);
    if (!strcmp(id, "S3")) return case_S3(emit);
    if (!strcmp(id, "T3")) return case_T3(emit);
    int r = find_r(id);
    if (r >= 0) return case_rprobe(r, emit);
    fprintf(stderr, "unknown case %s\n", id);
    exit(2);
}
static const char* all_cases[] = { "N3","S3","T3" };
#define NCASES 3

int main(int argc, char** argv) {
    if (argc >= 2 && !strcmp(argv[1], "list")) {
        for (int i = 0; i < NCASES; i++) puts(all_cases[i]);
        for (int i = 0; i < NR; i++) puts(iprobes[i].id);
        return 0;
    }
    if (argc >= 3 && !strcmp(argv[1], "bench")) {
        int n = argc >= 4 ? atoi(argv[3]) : 200000;
        for (int i = 0; i < n; i++) run_case(argv[2], 0);
        printf("sink=%d\n", (int)g_sink);
        return 0;
    }
    if (argc < 2) { fprintf(stderr, "usage: %s <case|bench|list>\n", argv[0]); return 2; }
    run_case(argv[1], 1);
    return 0;
}
