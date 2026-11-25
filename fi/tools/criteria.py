#!/usr/bin/env python3
"""Layered fault criteria for the mini-phys-lab injection suite.

Notation
--------
  r        return code of the injected call
  O        caller-provided pure out-params (Contact c, float loss):
           the caller reads them after *any* return, so they must be defined.
  S        authoritative in/out state (V3 v, Body b): a rejected call must be atomic.
  fin(x)   x is neither NaN nor +-Inf
  deg      the deliberately degenerate numeric input of the case
  ord/perm ordered and permuted runs of the same independent calls

Each criterion is a conjunction of checks ordered shallow -> deep. That order is
the point: a shallow check that passes on a faulty call is a *missed layer*.
A fault is "masked at layer L" when L passes and a deeper check is the first to
fail.

  C_num   = in_codes(r)        # the call claims a legitimate outcome (shallow)
            and finite(O)      # reported numbers are finite / in range
            and reject(deg)    # a degenerate number was actually refused

  C_state = defined(O/S)       # all reported fields finite, nothing NaN
            and cleared(O)     # on miss/reject no prior-frame value survives
            and atomic(S)      # a rejected in/out call left state bit-unchanged

  C_time  = ret_eq             # return codes agree across orders
            and fin_eq         # outputs finite in both orders
            and val_eq         # every observable bitwise equal across orders

A case PASSES iff its own-class criterion passes. Cross-class checks are run
only when the record actually carries that criterion's inputs; they surface
latent defects incidentally visible on the same call.
"""
import math

CODES = {"collide": {-2, -1, 0, 1}, "solve": {-1, 0, 1}, "integrate": {-1, 0, 1}}

CASE_META = {
    "N1": {"module": "collide",   "class": "num",   "inject": "sphere radius = 0"},
    "N2": {"module": "solve",     "class": "num",   "inject": "contact normal = (0,0,0)"},
    "N3": {"module": "integrate", "class": "num",   "inject": "dt = +Inf"},
    "S1": {"module": "collide",   "class": "state", "inject": "prior-frame Contact seed"},
    "S2": {"module": "solve",     "class": "state", "inject": "prior-frame loss seed"},
    "S3": {"module": "integrate", "class": "state", "inject": "prior-frame prev-position seed"},
    "T1": {"module": "collide",   "class": "time",  "inject": "4-way call permutation"},
    "T2": {"module": "solve",     "class": "time",  "inject": "3-way call permutation"},
    "T3": {"module": "integrate", "class": "time",  "inject": "3-way call permutation"},
}

# which records expose a given criterion's inputs
SIBLING_OK = {
    "num":   {"N1", "N2", "N3"},
    "state": {"N1", "S1", "S2", "S3"},
    "time":  {"T1", "T2", "T3"},
}


def num(v):
    if isinstance(v, str):
        return {"NaN": math.nan, "Inf": math.inf, "-Inf": -math.inf}.get(v, math.nan)
    return float(v)


def fin(v):
    try:
        return math.isfinite(num(v))
    except (TypeError, ValueError):
        return False


def finite_all(obj):
    if isinstance(obj, dict):
        return all(finite_all(x) for x in obj.values())
    if isinstance(obj, list):
        return all(finite_all(x) for x in obj)
    if isinstance(obj, str) and obj in ("NaN", "Inf", "-Inf"):
        return False
    return True


def eq(a, b, eps=1e-9):
    if isinstance(a, str) or isinstance(b, str):
        return a == b
    return abs(num(a) - num(b)) <= eps


def vzero(v, eps=1e-9):
    return all(abs(num(x)) <= eps for x in v)


def contact_cleared(c):
    return (c["hit"] == 0 and eq(c["depth"], 0.0) and vzero(c["n"]) and vzero(c["p"]))


class Verdict:
    def __init__(self, criterion):
        self.criterion = criterion
        self.checks = []

    def check(self, name, ok, detail=""):
        self.checks.append((name, bool(ok), detail))

    @property
    def ok(self):
        return all(ok for _, ok, _ in self.checks)

    def first_fail(self):
        for i, (n, ok, _) in enumerate(self.checks):
            if not ok:
                return i, n
        return None, None

    def all_failed(self):
        return [n for n, ok, _ in self.checks if not ok]

    def passed_above_failure(self):
        i, _ = self.first_fail()
        if i is None:
            return []
        return [n for n, ok, _ in self.checks[:i] if ok]

    def to_dict(self):
        return {"criterion": self.criterion, "pass": self.ok,
                "checks": [{"check": n, "pass": ok, "detail": d}
                           for n, ok, d in self.checks]}


# ------------------------------------------------------------------ C_num

def c_num(meta, rec):
    m = meta["module"]
    v = Verdict("C_num")
    r = rec["ret"]
    v.check("in_codes", r in CODES[m], f"ret={r}; documented={sorted(CODES[m])}")
    if m == "collide":
        c = rec["after"]
        v.check("finite_output",
                fin(c["depth"]) and all(fin(x) for x in c["n"]) and all(fin(x) for x in c["p"]),
                f"depth={c['depth']} n={c['n']} p={c['p']}")
        v.check("reject_deg", r < 0, f"radius=0 -> ret={r} (want <0)")
    elif m == "solve":
        loss = rec.get("loss")
        fin_ok = all(fin(x) for x in rec["v"]) and (loss == "null" or fin(loss))
        rng_ok = (loss == "null") or (num(loss) >= -1e-6)
        v.check("finite_output", fin_ok and rng_ok, f"v={rec['v']} loss={loss}")
        v.check("reject_deg", r < 0, f"|normal|=0 -> ret={r} (want <0)")
    else:
        b = rec["after"]
        v.check("finite_output",
                all(fin(b[f]) for f in ("x", "y", "px", "py", "vx", "vy")), f"body={b}")
        v.check("reject_deg", r < 0, f"dt=+Inf -> ret={r} (want <0)")
    return v


# ---------------------------------------------------------------- C_state

def c_state(meta, rec):
    m = meta["module"]
    v = Verdict("C_state")
    r = rec["ret"]
    if m == "collide":
        c = rec["after"]
        v.check("defined", finite_all(c), f"hit={c['hit']} depth={c['depth']} n={c['n']} p={c['p']}")
        cleared = contact_cleared(c)
        v.check("cleared", (r > 0) or cleared,
                "ok" if (r > 0 or cleared) else
                f"r={r} leaves prior contact in out-param: n={c['n']} p={c['p']}")
        v.check("atomic", True, "n/a: pure out-param")
    elif m == "solve":
        loss = rec["after"]
        v.check("defined", loss == "null" or fin(loss), f"loss={loss}")
        cleared = (loss != "null") and fin(loss) and abs(num(loss)) <= 1e-12
        v.check("cleared", (r >= 0) or cleared,
                "ok" if (r >= 0 or cleared) else
                f"r={r} leaves prior loss={loss} (before={rec['before']})")
        v.check("atomic", True, "n/a: pure out-param")
    else:
        b = rec["after"]
        v.check("defined", all(fin(b[f]) for f in ("x", "y", "px", "py", "vx", "vy")),
                "body fields finite" if all(fin(b[f]) for f in ("x", "y", "px", "py", "vx", "vy")) else f"body={b}")
        v.check("cleared", True, "n/a: authoritative state, not an out-param slot")
        if r < 0:
            same = all(eq(rec["before"][f], b[f]) for f in ("x", "y", "px", "py", "vx", "vy"))
            v.check("atomic", same, "rejected step unchanged" if same else "rejected step mutated body")
        else:
            v.check("atomic", True, "accepted step")
    return v


# ----------------------------------------------------------------- C_time

def c_time(meta, rec):
    v = Verdict("C_time")
    m = meta["module"]
    keys = {"collide": ["c"], "solve": ["v", "loss"], "integrate": ["body"]}[m]
    o, p = rec["ordered"], rec["permuted"]
    v.check("ret_eq", all(a.get("ret") == b.get("ret") for a, b in zip(o, p)),
            "return codes equal across orders")
    v.check("fin_eq", finite_all(o) and finite_all(p), "outputs finite in both orders")
    diffs = []
    for a, b in zip(o, p):
        for k in keys:
            if a[k] != b[k]:
                diffs.append((a.get("id"), k, a[k], b[k]))
    v.check("val_eq", not diffs and num(rec.get("maxdiff", 0)) == 0.0,
            "identical" if not diffs else f"divergence {diffs[:3]}")
    return v


CRIT = {"num": c_num, "state": c_state, "time": c_time}


def evaluate(case, rec):
    meta = CASE_META[case]
    own_name = meta["class"]
    own = CRIT[own_name](meta, rec)

    siblings = {}
    for name, fn in CRIT.items():
        if name == own_name:
            continue
        if case in SIBLING_OK[name]:
            siblings[name] = fn(meta, rec)

    cross = []
    for name, sib in siblings.items():
        if own.ok and not sib.ok:
            i, caught = sib.first_fail()
            cross.append({
                "caught_by": name, "caught_at_check": caught,
                "detail": sib.checks[i][2],
                "note": "own criterion passed; a latent defect is visible on the same call",
            })

    result = {"case": case, **meta, "own": own.to_dict(),
              "siblings": {k: s.to_dict() for k, s in siblings.items()},
              "cross_class": cross}

    if own.ok:
        result.update(verdict="PASS", caught_at=None, layers_that_missed=[],
                      failed_checks=[])
    else:
        i, caught = own.first_fail()
        result.update(verdict="FAIL", caught_at=caught,
                      layers_that_missed=own.passed_above_failure(),
                      failed_checks=own.all_failed(),
                      failure_detail=own.checks[i][2])
    return result
