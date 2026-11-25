#!/usr/bin/env python3
"""Run the 9-case injection suite and emit verdicts.

Each case runs as its own process through the plain (non-instrumented) adapter,
so injection and measurement are separated. Criteria are applied in
tools/criteria.py. Output: fi/results/verdicts.json + a human table on stdout.

Exit code: 0 if every case's own-class criterion passed, else 1.
"""
import json, os, subprocess, sys

FI = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BUILD = os.path.join(FI, "build")
RESULTS = os.path.join(FI, "results")
sys.path.insert(0, os.path.join(FI, "tools"))
import criteria  # noqa: E402

CORE = [("N1", "collide"), ("S1", "collide"), ("T1", "collide"),
        ("N2", "solve"),   ("S2", "solve"),   ("T2", "solve"),
        ("N3", "integrate"), ("S3", "integrate"), ("T3", "integrate")]


def run_case(module, case):
    exe = os.path.join(BUILD, module)
    out = subprocess.run([exe, case], check=True, capture_output=True, text=True).stdout.strip()
    return json.loads(out)


def main():
    os.makedirs(RESULTS, exist_ok=True)
    verdicts = []
    for case, module in CORE:
        rec = run_case(module, case)
        v = criteria.evaluate(case, rec)
        verdicts.append(v)

    with open(os.path.join(RESULTS, "verdicts.json"), "w") as fh:
        json.dump(verdicts, fh, indent=2)

    # human table
    print(f"{'case':4s} {'class':6s} {'module':10s} {'inject':30s} {'verdict':8s} "
          f"{'caught_at':18s} layers_that_missed")
    print("-" * 118)
    npass = 0
    for v in verdicts:
        if v["verdict"] == "PASS":
            npass += 1
        missed = ",".join(v.get("layers_that_missed", [])) or "-"
        caught = v.get("caught_at") or "-"
        print(f"{v['case']:4s} {v['class']:6s} {v['module']:10s} {v['inject']:30s} "
              f"{v['verdict']:8s} {caught:18s} {missed}")
    print("-" * 118)
    print(f"PASS {npass}/9  FAIL {9-npass}/9")

    # masking analysis
    # A masked fault = the case's own-class criterion passed, yet a defect is
    # visible on the same call under another criterion. We name the deep check
    # that fired AND the own-class check that was blind to it.
    cross = [v for v in verdicts if v.get("cross_class")]
    print(f"\ncross-class masking findings: {len(cross)}")
    for v in cross:
        for m in v["cross_class"]:
            print(f"  {v['case']} ({v['class']} fault): own criterion PASSED; "
                  f"caught by {m['caught_by']} at '{m['caught_at_check']}' -- {m['detail']}")

    # intra-criterion masking: for every FAIL, the shallow checks that passed
    # while the fault was already present are the layers that failed to catch it.
    print("\nin-criterion layer failures (guards that were blind before one fired):")
    for v in verdicts:
        if v["verdict"] == "FAIL":
            missed = ", ".join(v["layers_that_missed"]) or "(none: first check fired)"
            print(f"  {v['case']}: caught at '{v['caught_at']}'; "
                  f"blind layers above = [{missed}]")

    return 0 if npass == 9 else 1


if __name__ == "__main__":
    sys.exit(main())
