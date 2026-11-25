#!/usr/bin/env python3
"""Coverage + timing collector for fault-injection adapters.

A "path" here is one coverage unit: an executable source line or one arm of a
conditional branch, attributed to the function it lives in. That is the standard
reachable-unit proxy for path coverage; we state it explicitly rather than claim
true path enumeration.

Scope: a module's functional functions (clampf/closest/test, resolve, step).
The module's own self-test (check + the embedded main) is excluded: no injection
point drives it, so counting it would understate coverage for everyone.

Usage:
  covtool.py matrix <module> <case>...
      -> JSON {case: {"lines": [...], "branches": [[l,a],...],
                      "funcs": {fn: {"lines": [...], "branches": [...]}}}}
  covtool.py declared <module>
      -> JSON {fn: {"lines": [...], "branches": [...]}} for every function
  covtool.py time <module> <case> [n1] [n2]
      -> JSON {"case","percall_ns","n1","n2"}
"""
import os, re, subprocess, sys, json, tempfile, shutil, glob

FI = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # .../fi
BUILD = os.path.join(FI, "build")

LINE_RE = re.compile(r'^\s*([^:]+):\s*(\d+):(.*)$')
FUNC_RE = re.compile(r'^function (\S+) called')
BRANCH_RE = re.compile(r'^branch\s+(\d+)\s+(?:taken\s+(\d+)|never executed)')

SELF_TEST_FUNCS = {"check"}
def is_selftest(fn):
    return fn in SELF_TEST_FUNCS or fn.endswith("module_main")


def parse_gcov(path):
    """Return {func: {"lines": set, "branches": set[(line,arm)]}} for one .gcov."""
    data = {}
    cur_fn, cur_line = None, None
    with open(path, "r", errors="replace") as fh:
        for raw in fh:
            fm = FUNC_RE.match(raw)
            if fm:
                cur_fn = fm.group(1)
                data.setdefault(cur_fn, {"lines": set(), "branches": set()})
                continue
            m = LINE_RE.match(raw.rstrip("\n"))
            if m:
                cnt, lineno = m.group(1).strip(), int(m.group(2))
                cur_line = lineno
                if lineno and cnt not in ("-", "") and cur_fn:
                    if cnt != "#####":
                        try:
                            if int(cnt.rstrip("*^$")) > 0:
                                data[cur_fn]["lines"].add(lineno)
                        except ValueError:
                            pass
                continue
            bm = BRANCH_RE.match(raw.strip())
            if bm and cur_line is not None and cur_fn:
                arm = int(bm.group(1))
                if bm.group(2) is not None and int(bm.group(2)) > 0:
                    data[cur_fn]["branches"].add((cur_line, arm))
    return data


def _gcov_for(module, case):
    exe = os.path.join(BUILD, f"{module}_cov")
    tmp = tempfile.mkdtemp(prefix=f"cov_{module}_{case}_")
    env = dict(os.environ, GCOV_PREFIX=tmp, GCOV_PREFIX_STRIP="0")
    subprocess.run([exe, case], check=True, env=env,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    gcda = glob.glob(os.path.join(tmp, "**", "*.gcda"), recursive=True)
    if not gcda:
        shutil.rmtree(tmp, ignore_errors=True)
        raise RuntimeError(f"no .gcda produced for {module}/{case}")
    gcno = os.path.join(BUILD, os.path.basename(gcda[0])[:-5] + ".gcno")
    if os.path.exists(gcno):
        shutil.copy(gcno, os.path.join(os.path.dirname(gcda[0]), os.path.basename(gcno)))
    subprocess.run(["gcov", "-b", "-c", gcda[0]], cwd=tmp, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    src = os.path.join(tmp, f"{module}.c.gcov")
    if not os.path.exists(src):
        shutil.rmtree(tmp, ignore_errors=True)
        raise RuntimeError(f"no module .gcov for {module}/{case}")
    return tmp, src


def coverage(module, case):
    tmp, src = _gcov_for(module, case)
    try:
        return parse_gcov(src)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def list_cases(module):
    exe = os.path.join(BUILD, module)
    out = subprocess.run([exe, "list"], check=True, capture_output=True, text=True).stdout
    return out.split()


def _flatten(data, include_selftest=False):
    lines, branches, funcs = set(), set(), {}
    for fn, d in data.items():
        if not include_selftest and is_selftest(fn):
            continue
        funcs[fn] = {"lines": sorted(d["lines"]), "branches": sorted(d["branches"])}
        lines |= d["lines"]
        branches |= d["branches"]
    return lines, branches, funcs


def declared(module):
    """Declared (whether taken or not) units: run one case, parse gcov source once."""
    exe = os.path.join(BUILD, f"{module}_cov")
    tmp = tempfile.mkdtemp(prefix=f"decl_{module}_")
    env = dict(os.environ, GCOV_PREFIX=tmp, GCOV_PREFIX_STRIP="0")
    subprocess.run([exe, list_cases(module)[0]], check=True, env=env,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    gcda = glob.glob(os.path.join(tmp, "**", "*.gcda"), recursive=True)[0]
    gcno = os.path.join(BUILD, os.path.basename(gcda)[:-5] + ".gcno")
    shutil.copy(gcno, os.path.join(os.path.dirname(gcda), os.path.basename(gcno)))
    subprocess.run(["gcov", "-b", "-c", gcda], cwd=tmp, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    src = os.path.join(tmp, f"{module}.c.gcov")
    # re-parse keeping only the *declared* set: executable lines + all branch arms
    data = {}
    cur_fn, cur_line = None, None
    for raw in open(src, errors="replace"):
        fm = FUNC_RE.match(raw)
        if fm:
            cur_fn = fm.group(1); data.setdefault(cur_fn, {"lines": set(), "branches": set()}); continue
        m = LINE_RE.match(raw.rstrip("\n"))
        if m:
            cnt, lineno = m.group(1).strip(), int(m.group(2)); cur_line = lineno
            if lineno and cnt not in ("-", "") and cur_fn:
                data[cur_fn]["lines"].add(lineno)   # declared regardless of count
            continue
        bm = BRANCH_RE.match(raw.strip())
        if bm and cur_line is not None and cur_fn:
            data[cur_fn]["branches"].add((cur_line, int(bm.group(1))))
    shutil.rmtree(tmp, ignore_errors=True)
    return data


def main():
    if len(sys.argv) < 3:
        print(__doc__); return 2
    cmd, module = sys.argv[1], sys.argv[2]

    if cmd == "matrix":
        out = {}
        for c in sys.argv[3:]:
            lines, branches, funcs = _flatten(coverage(module, c))
            out[c] = {"lines": sorted(lines), "branches": sorted(branches), "funcs": funcs}
        json.dump(out, sys.stdout, indent=1); print(); return 0

    if cmd == "declared":
        d = declared(module)
        lines, branches, funcs = _flatten(d)
        json.dump({"funcs": funcs, "lines": sorted(lines), "branches": sorted(branches),
                   "n_units": len(lines) + len(branches)}, sys.stdout, indent=1)
        print(); return 0

    if cmd == "time":
        import time
        case = sys.argv[3]
        n1 = int(sys.argv[4]) if len(sys.argv) > 4 else 100000
        n2 = int(sys.argv[5]) if len(sys.argv) > 5 else 200000
        exe = os.path.join(BUILD, module)

        def run(n):
            best = float("inf")
            for _ in range(5):
                t0 = time.perf_counter()
                subprocess.run([exe, "bench", case, str(n)], check=True,
                               stdout=subprocess.DEVNULL)
                best = min(best, time.perf_counter() - t0)
            return best
        t1, t2 = run(n1), run(n2)
        percall = (t2 - t1) / (n2 - n1)
        json.dump({"module": module, "case": case, "percall_ns": percall * 1e9,
                   "n1": n1, "n2": n2}, sys.stdout); print(); return 0

    print(f"unknown command {cmd}", file=sys.stderr); return 2


if __name__ == "__main__":
    sys.exit(main())
