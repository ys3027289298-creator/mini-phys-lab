#!/usr/bin/env python3
"""Map coverage units (source lines) to named code paths.

A "path" here is a semantically named decision region inside a functional
function, e.g. collide's radius-reject guard, the surface-normal branch, or the
axis-selection fan inside the box. The mapping is by source line, which is
stable for this commit. Branch arms remain the precise second-level signal in
coverage_matrix.json; this lexicon is the human-facing grouping.

Covered line -> path label; a point "covers" a path iff it covers any of its
lines. Paths are library-side only (self-tests excluded), matching covtool.
"""
# (module, line) -> path label
LEXICON = {
    "collide": {
        13: "clamp:entry",
        14: "clamp:low", 15: "clamp:high", 16: "clamp:mid",
        20: "closest:project", 22: "closest:project", 23: "closest:project",
        24: "closest:project", 25: "closest:project",
        28: "entry:test",
        29: "guard:radius", 30: "guard:box",
        31: "compute:delta", 32: "compute:delta", 33: "compute:delta",
        34: "reset:out", 35: "reset:out",
        36: "decision:miss",
        37: "decision:surface", 38: "surface:distance", 39: "surface:normal",
        40: "surface:depth",
        43: "inside:shallow", 44: "inside:shallow", 45: "inside:shallow", 46: "inside:shallow",
        47: "inside:axisY", 48: "inside:axisZ",
        49: "inside:normal-reset",
        50: "inside:signX", 51: "inside:signY", 52: "inside:signZ",
        53: "inside:depth",
        55: "emit:point", 56: "emit:hit", 57: "return:hit",
    },
    "solve": {
        15: "entry:resolve",
        16: "guard:restitution", 17: "guard:friction",
        19: "compute:vn", 20: "decision:separating",
        21: "separating:zero-loss", 22: "separating:return",
        25: "compute:tangent", 26: "compute:tangent",
        29: "normal:new-vn",
        31: "friction:reduce", 32: "friction:clip", 33: "friction:scale",
        35: "apply:velocity", 36: "apply:velocity", 37: "apply:velocity",
        39: "out:loss", 40: "return:resolved",
    },
    "integrate": {
        9: "entry:step",
        10: "guard:dt-radius",
        11: "save:prev",
        12: "integrate:velocity", 13: "integrate:velocity",
        14: "integrate:position", 15: "integrate:position",
        17: "compute:penetration", 18: "decision:penetration",
        19: "correct:position", 20: "correct:stop-velocity",
        22: "return:status",
    },
}


def paths_for_lines(module, lines):
    lex = LEXICON[module]
    return sorted({lex[l] for l in lines if l in lex})


def all_paths(module, declared_lines):
    return paths_for_lines(module, declared_lines)
