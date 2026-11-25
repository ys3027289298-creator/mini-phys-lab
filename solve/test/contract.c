/* contract test: the type argument is mandatory and unforgeable.

   Run via `make contract`. Each case is compiled in isolation; the verdict
   is decided by whether the compiler emits an error. A case that should be
   rejected and compiles clean is a FAILURE - it means the design stopped
   protecting against the incident. Kept out of the normal build because two
   of the cases are deliberate compile errors. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define INCDIR_SOLVE   (argc > 1 ? argv[1] : "src")
#define INCDIR_INTEG   (argc > 2 ? argv[2] : "../integrate/src")

typedef struct { const char* name; const char* src; int want_error; } Case;

int main(int argc, char** argv) {
    static const char* footer =
        "int main(void){return 0;}\n";

    /* -I paths: the two source dirs are fixed relative to the repo root; the
       runner is invoked from solve/ so src/ and ../integrate/src/ resolve. */
    const char* inc_solve = (argc > 1 ? argv[1] : "src");
    const char* inc_integ = (argc > 2 ? argv[2] : "../integrate/src");

    static const Case cases[] = {
        /* --- must be rejected: type omitted (the incident shape) --- */
        { "omit_type_integrate", "#include \"integrate.h\"\n"
              "int f(Body* b){ return step(b, 0.1f, 0, -9.8f, 0, 0.5f); }\n", 1 },
        { "omit_type_solve",     "#include \"solve.h\"\n"
              "int f(V3* v, V3 n, const Material* m){ float l; return resolve(v, n, m, &l); }\n", 1 },
        /* --- must be rejected: literal standing in for the kind --- */
        { "literal_zero_integrate", "#include \"integrate.h\"\n"
              "int f(Body* b){ return step(0, b, 0.1f, 0, -9.8f, 0, 0.5f); }\n", 1 },
        { "literal_one_solve",      "#include \"solve.h\"\n"
              "int f(V3* v, V3 n, const Material* m){ float l; return resolve(1, v, n, m, &l); }\n", 1 },
        /* --- must compile with zero diagnostics: both kinds, both modules --- */
        { "correct_dynamic_solve",  "#include \"solve.h\"\n"
              "int f(V3* v, V3 n, const Material* m){ float l; return resolve(BODY_DYNAMIC, v, n, m, &l); }\n", 0 },
        { "correct_kinematic_solve","#include \"solve.h\"\n"
              "int f(V3* v, V3 n, const Material* m){ float l; return resolve(BODY_KINEMATIC, v, n, m, &l); }\n", 0 },
        { "correct_kinematic_step", "#include \"integrate.h\"\n"
              "int f(Body* b){ return step(BODY_KINEMATIC, b, 0.1f, 0, -9.8f, 0, 0.5f); }\n", 0 },
        /* --- both headers in one TU (a scene includes both) --- */
        { "both_headers_one_tu",    "#include \"solve.h\"\n#include \"integrate.h\"\n"
              "int f(Body* b){ return step(BODY_DYNAMIC, b, 0.1f, 0, -9.8f, 0, 0.5f); }\n", 0 },
    };
    const int n = (int)(sizeof cases / sizeof cases[0]);

    int failures = 0;
    char path[256], cmd[1024];
    for (int i = 0; i < n; i++) {
        const Case* c = &cases[i];
        snprintf(path, sizeof path, "/tmp/contract_%s.c", c->name);
        FILE* fp = fopen(path, "w");
        if (!fp) { fprintf(stderr, "cannot write %s\n", path); return 2; }
        fputs(c->src, fp); fputs(footer, fp); fclose(fp);

        /* compile with the project's warning set; -Werror turns every warning
           into an error, so "clean" means truly silence-free. */
        snprintf(cmd, sizeof cmd,
                 "%s -std=c11 -Wall -Wextra -Werror -Wswitch-enum -Wpedantic "
                 "-I %s -I %s -c %s -o /dev/null 2>/dev/null",
                 (getenv("CC") ? getenv("CC") : "cc"), inc_solve, inc_integ, path);
        int rc = system(cmd);
        int got_error = (rc != 0);
        int ok = got_error == c->want_error;
        printf("%s=%d\n", c->name, ok ? 1 : 0);
        if (!ok) failures++;
        remove(path);
    }
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
