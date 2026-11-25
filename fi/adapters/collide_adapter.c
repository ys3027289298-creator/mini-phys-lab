/* Fault-injection adapter for collide.
 *
 * Includes the unmodified module so its static internals (test/closest/clampf)
 * are reachable, then drives them with injected faults. Each case runs in its
 * own process so gcov attribution stays isolated. `emit=0` is the timing path:
 * same work, no serialisation, result folded into g_sink.
 */
#define main collide_module_main
#include "../../collide/src/collide.c"
#undef main
#include <stdlib.h>

static volatile int g_sink;

static void jf(float v) {
    if (isnan(v)) fputs("\"NaN\"", stdout);
    else if (isinf(v)) fputs(v > 0 ? "\"Inf\"" : "\"-Inf\"", stdout);
    else printf("%.9g", v);
}
static void jv(V3 v) { putchar('['); jf(v.x); putchar(','); jf(v.y); putchar(','); jf(v.z); putchar(']'); }
static void jcontact(Contact* c) {
    printf("{\"hit\":%d,\"depth\":", c->hit); jf(c->depth);
    printf(",\"n\":"); jv(c->normal); printf(",\"p\":"); jv(c->point); printf("}");
}

static const Aabb kbox = { {0,0,0}, {1,1,1} };

/* ---- NUM: zero radius -------------------------------------------------- */
static void case_N1(int emit) {
    Sphere s = { {1.4f,0.5f,0.5f}, 0.0f };
    Contact c = { {0.7f,-0.7f,0.1f}, 9.0f, {7,8,9}, 1 };   /* prior-frame sentinel */
    Contact before = c;
    int ret = test(&s, &kbox, &c);
    g_sink += ret;
    if (!emit) return;
    printf("{\"module\":\"collide\",\"case\":\"N1\",\"ret\":%d,\"before\":", ret);
    jcontact(&before); printf(",\"after\":"); jcontact(&c); printf("}\n");
}

/* ---- STATE: residual contact then a miss ------------------------------- */
static void case_S1(int emit) {
    Sphere far = { {5,5,5}, 0.5f };
    Contact c = { {0.7f,-0.7f,0.1f}, 0.123f, {7,8,9}, 1 }; /* prior-frame contact */
    Contact before = c;
    int ret = test(&far, &kbox, &c);
    g_sink += ret;
    if (!emit) return;
    printf("{\"module\":\"collide\",\"case\":\"S1\",\"ret\":%d,\"before\":", ret);
    jcontact(&before); printf(",\"after\":");
    jcontact(&c); printf("}\n");
}

/* ---- TIMING: order permutation over independent spheres ---------------- */
static const Sphere t1_spheres[4] = {
    { {1.4f,0.5f,0.5f}, 0.5f },   /* face   */
    { {1.2f,1.2f,1.2f}, 0.5f },   /* corner */
    { {0.5f,0.5f,0.5f}, 0.25f },  /* inside */
    { {1.5f,0.5f,0.5f}, 0.5f },   /* grazing */
};
static void t1_run(const int* order, int out_ret[4], Contact out_c[4]) {
    for (int k = 0; k < 4; k++) {
        int id = order[k];
        Contact c = { {0,0,0}, 0.0f, {0,0,0}, 0 };
        out_ret[id] = test(&t1_spheres[id], &kbox, &c);
        out_c[id] = c;
    }
}
static void case_T1(int emit) {
    static const int ord[4]  = {0,1,2,3};
    static const int perm[4] = {3,2,1,0};
    int oret[4], pret[4]; Contact oc[4], pc[4];
    t1_run(ord, oret, oc);
    t1_run(perm, pret, pc);
    float maxdiff = 0.0f;
    for (int i = 0; i < 4; i++) {
        float d = fabsf(oc[i].depth - pc[i].depth);
        d = fmaxf(d, fabsf(oc[i].normal.x - pc[i].normal.x));
        d = fmaxf(d, fabsf(oc[i].normal.y - pc[i].normal.y));
        d = fmaxf(d, fabsf(oc[i].normal.z - pc[i].normal.z));
        if (oret[i] != pret[i] || oc[i].hit != pc[i].hit) d = INFINITY;
        maxdiff = fmaxf(maxdiff, d);
        g_sink += oret[i] + pret[i];
    }
    if (!emit) return;
    printf("{\"module\":\"collide\",\"case\":\"T1\",\"maxdiff\":");
    jf(maxdiff);
    printf(",\"ordered\":[");
    for (int i = 0; i < 4; i++) { if (i) putchar(','); printf("{\"id\":%d,\"ret\":%d,\"c\":", i, oret[i]); jcontact(&oc[i]); putchar('}'); }
    printf("],\"permuted\":[");
    for (int i = 0; i < 4; i++) { if (i) putchar(','); printf("{\"id\":%d,\"ret\":%d,\"c\":", i, pret[i]); jcontact(&pc[i]); putchar('}'); }
    printf("]}\n");
}

/* ---- extra candidate injection points (coverage experiment) ------------ */
typedef struct { const char* id; V3 c; float r; Aabb b; } Probe;
#define KBOX ((Aabb){ {0,0,0}, {1,1,1} })
static const Probe probes[] = {
    {"P_face",      {1.4f,0.5f,0.5f}, 0.5f, KBOX},
    {"P_graz",      {1.5f,0.5f,0.5f}, 0.5f, KBOX},
    {"P_corner",    {1.2f,1.2f,1.2f}, 0.5f, KBOX},
    {"P_low",       {-0.5f,0.5f,0.5f},0.5f, KBOX},   /* centre below min: clamp low */
    {"P_insideXl",  {0.1f,0.5f,0.5f}, 0.25f,KBOX},   /* centre inside, shallow x, below mid */
    {"P_insideXh",  {0.9f,0.5f,0.5f}, 0.25f,KBOX},   /* centre inside, shallow x, above mid */
    {"P_insideYl",  {0.5f,0.1f,0.5f}, 0.25f,KBOX},
    {"P_insideYh",  {0.5f,0.9f,0.5f}, 0.25f,KBOX},
    {"P_insideZl",  {0.5f,0.5f,0.1f}, 0.25f,KBOX},
    {"P_insideZh",  {0.5f,0.5f,0.9f}, 0.25f,KBOX},
    {"P_miss",      {5,5,5},          0.5f, KBOX},
    {"P_missLow",   {-5,-5,-5},       0.5f, KBOX},
    {"P_inf",       {1.4f,0.5f,0.5f}, INFINITY, KBOX},
    {"P_nanr",      {1.4f,0.5f,0.5f}, 0.0f/0.0f, KBOX},
    {"P_invAll",    {0.5f,0.5f,0.5f}, 0.5f, { {1,1,1}, {0,0,0} }},
    {"P_invY",      {0.5f,0.5f,0.5f}, 0.5f, { {0,1,0}, {1,0,1} }},
    {"P_invZ",      {0.5f,0.5f,0.5f}, 0.5f, { {0,0,1}, {1,1,0} }},
};
#define NPRB ((int)(sizeof(probes)/sizeof(probes[0])))

static void case_probe(int idx, int emit) {
    const Probe* p = &probes[idx];
    Contact c = { {0,0,0}, 0.0f, {0,0,0}, 0 };
    int ret = test(&(Sphere){p->c, p->r}, &p->b, &c);
    g_sink += ret;
    if (!emit) return;
    printf("{\"module\":\"collide\",\"case\":\"%s\",\"ret\":%d,\"after\":", p->id, ret);
    jcontact(&c); printf("}\n");
}

/* ---- supplementary diagnostic: order dependence via a SHARED slot ------
 * The 9-case timing tests use independent slots, where these pure functions
 * cannot be order-dependent. This diagnostic reuses one Contact across
 * [hit, miss] in both orders to expose residual-induced order dependence.
 * It is evidence for the masking analysis, not one of the 9 mandated cases. */
static void case_X_alias(int emit) {
    Sphere face = { {1.4f,0.5f,0.5f}, 0.5f };
    Sphere miss = { {5,5,5}, 0.5f };
    Contact a = { {0,0,0}, 0.0f, {0,0,0}, 0 };
    Contact b = { {0,0,0}, 0.0f, {0,0,0}, 0 };
    test(&face, &kbox, &a); test(&miss, &kbox, &a);   /* ordered:  hit, miss   */
    test(&miss, &kbox, &b); test(&face, &kbox, &b);   /* permuted: miss, hit   */
    g_sink += a.hit + b.hit;
    if (!emit) return;
    printf("{\"module\":\"collide\",\"case\":\"X_alias\",\"ordered\":"); jcontact(&a);
    printf(",\"permuted\":"); jcontact(&b);
    printf(",\"equal\":%d}\n", a.hit == b.hit && a.depth == b.depth &&
           a.normal.x == b.normal.x && a.normal.y == b.normal.y &&
           a.normal.z == b.normal.z && a.point.x == b.point.x &&
           a.point.y == b.point.y && a.point.z == b.point.z);
}

/* ---- dispatch ---------------------------------------------------------- */
static int find_probe(const char* id) {
    for (int i = 0; i < NPRB; i++) if (!strcmp(probes[i].id, id)) return i;
    return -1;
}
static void run_case(const char* id, int emit) {
    if (!strcmp(id, "N1")) return case_N1(emit);
    if (!strcmp(id, "S1")) return case_S1(emit);
    if (!strcmp(id, "T1")) return case_T1(emit);
    if (!strcmp(id, "X_alias")) return case_X_alias(emit);
    int p = find_probe(id);
    if (p >= 0) return case_probe(p, emit);
    fprintf(stderr, "unknown case %s\n", id);
    exit(2);
}
static const char* all_cases[] = { "N1","S1","T1" };
#define NCASES 3

int main(int argc, char** argv) {
    if (argc >= 2 && !strcmp(argv[1], "list")) {
        for (int i = 0; i < NCASES; i++) puts(all_cases[i]);
        for (int i = 0; i < NPRB; i++) puts(probes[i].id);
        return 0;
    }
    if (argc >= 3 && !strcmp(argv[1], "bench")) {
        const char* id = argv[2];
        int n = argc >= 4 ? atoi(argv[3]) : 200000;
        for (int i = 0; i < n; i++) run_case(id, 0);
        printf("sink=%d\n", (int)g_sink);
        return 0;
    }
    if (argc < 2) { fprintf(stderr, "usage: %s <case|bench|list>\n", argv[0]); return 2; }
    run_case(argv[1], 1);
    return 0;
}
