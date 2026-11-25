/* Fault-injection adapter for solve. Includes the unmodified module so the
 * static resolve() is reachable. One case per process for clean gcov. */
#define main solve_module_main
#include "../../solve/src/solve.c"
#undef main
#include <stdlib.h>

static volatile int g_sink;

static void jf(float v) {
    if (isnan(v)) fputs("\"NaN\"", stdout);
    else if (isinf(v)) fputs(v > 0 ? "\"Inf\"" : "\"-Inf\"", stdout);
    else printf("%.9g", v);
}
static void jv(V3 v) { putchar('['); jf(v.x); putchar(','); jf(v.y); putchar(','); jf(v.z); putchar(']'); }

/* ---- NUM: zero-length normal (degenerate contact normal) --------------- */
static void case_N2(int emit) {
    V3 n = { 0.0f, 0.0f, 0.0f };
    V3 v = { 0.0f, -10.0f, 0.0f };
    Material m = { 0.5f, 0.2f };
    float loss = -1.0f;
    int ret = resolve(&v, n, &m, &loss);
    float nlen = sqrtf(n.x*n.x + n.y*n.y + n.z*n.z);
    g_sink += ret;
    if (!emit) return;
    printf("{\"module\":\"solve\",\"case\":\"N2\",\"ret\":%d,\"normal_len\":", ret); jf(nlen);
    printf(",\"v\":"); jv(v); printf(",\"loss\":"); jf(loss); printf("}\n");
}

/* ---- STATE: residual out-param survives a rejected call ---------------- */
static void case_S2(int emit) {
    V3 n = { 0, 1, 0 };
    V3 v = { 0.0f, -1.0f, 0.0f };
    Material bad = { 1.5f, 0.2f };            /* restitution out of range -> -1 */
    float loss = 0.123f;                       /* prior-frame residual */
    float before = loss;
    int ret = resolve(&v, n, &bad, &loss);
    g_sink += ret;
    if (!emit) return;
    printf("{\"module\":\"solve\",\"case\":\"S2\",\"ret\":%d,\"before\":", ret); jf(before);
    printf(",\"after\":"); jf(loss); printf(",\"v\":"); jv(v); printf("}\n");
}

/* ---- TIMING: permute independent resolves ------------------------------ */
static const V3 t2_v[3] = { {0,-10,0}, {2,-2,0}, {3,0,1} };
static const Material t2_m[3] = { {0.5f,0.2f}, {0.0f,0.0f}, {0.8f,0.5f} };
static void t2_run(const int* order, int out_ret[3], V3 out_v[3], float out_loss[3]) {
    V3 n = { 0, 1, 0 };
    for (int k = 0; k < 3; k++) {
        int id = order[k];
        V3 v = t2_v[id];
        float loss = 0.0f;
        out_ret[id] = resolve(&v, n, &t2_m[id], &loss);
        out_v[id] = v;
        out_loss[id] = loss;
    }
}
static void case_T2(int emit) {
    static const int ord[3]  = {0,1,2};
    static const int perm[3] = {2,0,1};
    int oret[3], pret[3]; V3 ov[3], pv[3]; float ol[3], pl[3];
    t2_run(ord, oret, ov, ol);
    t2_run(perm, pret, pv, pl);
    float maxdiff = 0.0f;
    for (int i = 0; i < 3; i++) {
        float d = fabsf(ol[i] - pl[i]);
        d = fmaxf(d, fabsf(ov[i].x - pv[i].x));
        d = fmaxf(d, fabsf(ov[i].y - pv[i].y));
        d = fmaxf(d, fabsf(ov[i].z - pv[i].z));
        if (oret[i] != pret[i]) d = INFINITY;
        maxdiff = fmaxf(maxdiff, d);
        g_sink += oret[i] + pret[i];
    }
    if (!emit) return;
    printf("{\"module\":\"solve\",\"case\":\"T2\",\"maxdiff\":"); jf(maxdiff);
    printf(",\"ordered\":[");
    for (int i = 0; i < 3; i++) { if (i) putchar(','); printf("{\"id\":%d,\"ret\":%d,\"v\":", i, oret[i]); jv(ov[i]); printf(",\"loss\":"); jf(ol[i]); putchar('}'); }
    printf("],\"permuted\":[");
    for (int i = 0; i < 3; i++) { if (i) putchar(','); printf("{\"id\":%d,\"ret\":%d,\"v\":", i, pret[i]); jv(pv[i]); printf(",\"loss\":"); jf(pl[i]); putchar('}'); }
    printf("]}\n");
}

/* ---- extra candidate injection points (coverage experiment) ------------ */
typedef struct { const char* id; V3 v; V3 n; Material m; int null_out; } SProbe;
static const SProbe sprobes[] = {
    {"Q_bounce",   {0,-10,0},   {0,1,0}, {0.5f,0.2f}, 0},
    {"Q_sep",      {0,3,0},     {0,1,0}, {0.5f,0.2f}, 0},
    {"Q_fricclip", {4,-1,0},    {0,1,0}, {0.0f,1.0f}, 0},
    {"Q_fricpart", {4,-10,0},   {0,1,0}, {0.5f,0.2f}, 0},
    {"Q_dull",     {2,-2,0},    {0,1,0}, {0.0f,0.0f}, 0},
    {"Q_badrestHi",{0,-1,0},    {0,1,0}, {1.5f,0.2f}, 0},
    {"Q_badrestLo",{0,-1,0},    {0,1,0}, {-0.1f,0.2f}, 0},
    {"Q_badfricLo",{0,-1,0},    {0,1,0}, {0.5f,-0.1f}, 0},
    {"Q_badfricHi",{0,-1,0},    {0,1,0}, {0.5f,1.5f}, 0},
    {"Q_zeron",    {0,-10,0},   {0,0,0}, {0.5f,0.2f}, 0},
    {"Q_nanvel",   {0.0f/0.0f,-1,0}, {0,1,0}, {0.5f,0.2f}, 0},
    {"Q_vnzero",   {4,0,0},     {0,1,0}, {0.5f,0.2f}, 0},
    {"Q_tinytan",  {1e-12f,-10,0}, {0,1,0}, {0.5f,0.2f}, 0},
    {"Q_sepNull",  {0,3,0},     {0,1,0}, {0.5f,0.2f}, 1},   /* line 21, out==NULL */
    {"Q_hitNull",  {0,-10,0},   {0,1,0}, {0.5f,0.2f}, 1},   /* line 39, out==NULL */
};
#define NQ ((int)(sizeof(sprobes)/sizeof(sprobes[0])))

static void case_qprobe(int idx, int emit) {
    const SProbe* p = &sprobes[idx];
    V3 v = p->v; float loss = 0.0f;
    int ret = resolve(&v, p->n, &p->m, p->null_out ? NULL : &loss);
    g_sink += ret;
    if (!emit) return;
    printf("{\"module\":\"solve\",\"case\":\"%s\",\"ret\":%d,\"v\":", p->id, ret); jv(v);
    printf(",\"loss\":"); if (p->null_out) fputs("\"null\"", stdout); else jf(loss); printf("}\n");
}

static int find_q(const char* id) {
    for (int i = 0; i < NQ; i++) if (!strcmp(sprobes[i].id, id)) return i;
    return -1;
}
static void run_case(const char* id, int emit) {
    if (!strcmp(id, "N2")) return case_N2(emit);
    if (!strcmp(id, "S2")) return case_S2(emit);
    if (!strcmp(id, "T2")) return case_T2(emit);
    int q = find_q(id);
    if (q >= 0) return case_qprobe(q, emit);
    fprintf(stderr, "unknown case %s\n", id);
    exit(2);
}
static const char* all_cases[] = { "N2","S2","T2" };
#define NCASES 3

int main(int argc, char** argv) {
    if (argc >= 2 && !strcmp(argv[1], "list")) {
        for (int i = 0; i < NCASES; i++) puts(all_cases[i]);
        for (int i = 0; i < NQ; i++) puts(sprobes[i].id);
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
