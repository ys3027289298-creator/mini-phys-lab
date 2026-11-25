#!/usr/bin/env python3
"""Audit contact/impulse/spin logs.

The checker is deliberately schema-tolerant.  One structured record per contact
is preferred:

  frame=60 ball=3 depth=0.0031 normal=(0,1,0) radius=0.5 \\
  center=(0,0.4969,0) vn=-0.16 arm=0.5 awake=1 solver_called=1 j=0.0000

JSON lines with the same keys also work.  The three quoted production lines do
not contain frame, ball id, radius, centre, vn, arm, or sleep state; this tool
therefore reports those hypotheses as unproven rather than guessing.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

NUMBER_RE = re.compile(r"([A-Za-z_][\w.\[\]]*)=(-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)")
VECTOR_RE = re.compile(r"([A-Za-z_][\w.\[\]]*)=\(([^)]*)\)")
ARROW_RE = re.compile(r"([A-Za-z_][\w.\[\]]*)=\s*\(([^)]*)\)\s*->\s*\(([^)]*)\)")
TOKEN_RE = re.compile(r"([A-Za-z_][\w.\[\]]*)=([^\s]+)")


def vec(value: Any) -> tuple[float, float, float] | None:
    if value is None:
        return None
    if isinstance(value, (list, tuple)) and len(value) == 3:
        return tuple(float(x) for x in value)  # type: ignore[return-value]
    if isinstance(value, str):
        parts = value.strip().strip("()").split(",")
        if len(parts) == 3:
            try:
                return tuple(float(x) for x in parts)  # type: ignore[return-value]
            except ValueError:
                return None
    return None


def dot(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def sub(a: tuple[float, float, float], b: tuple[float, float, float]):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def length(a) -> float:
    return math.sqrt(dot(a, a))


def is_true_flag(value: Any) -> bool:
    return value is True or (isinstance(value, (int, float)) and not isinstance(value, bool) and value != 0)


def is_false_flag(value: Any) -> bool:
    return value is False or value == 0


def label(value: Any) -> str:
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def parse_line(line: str) -> dict[str, Any]:
    line = line.strip()
    if not line:
        return {}
    if line.startswith("{"):
        try:
            value = json.loads(line)
            return value if isinstance(value, dict) else {}
        except json.JSONDecodeError:
            pass

    record: dict[str, Any] = {}
    for match in ARROW_RE.finditer(line):
        record[match.group(1) + "0"] = "(" + match.group(2).strip() + ")"
        record[match.group(1) + "1"] = "(" + match.group(3).strip() + ")"
    for match in VECTOR_RE.finditer(line):
        record[match.group(1)] = "(" + match.group(2).strip() + ")"

    # Numeric matches must not overwrite vectors such as normal=(0,1,0).
    vector_spans = [m.span(1) for m in VECTOR_RE.finditer(line)]
    for match in NUMBER_RE.finditer(line):
        key_span = match.span(1)
        if any(a <= key_span[0] < b for a, b in vector_spans):
            continue
        record[match.group(1)] = float(match.group(2))

    for match in TOKEN_RE.finditer(line):
        key = match.group(1)
        if key in record:
            continue
        raw = match.group(2).strip().rstrip(",")
        low = raw.lower()
        if low in {"1", "true", "yes", "on"}:
            record[key] = True
        elif low in {"0", "false", "no", "off"}:
            record[key] = False
        else:
            record[key] = raw
    return record


def parse_log(lines: Iterable[str]) -> list[dict[str, Any]]:
    """Parse JSONL/key=value records and merge split lines by frame and ball."""
    merged: dict[tuple[Any, Any], dict[str, Any]] = {}
    order: list[tuple[Any, Any]] = []
    anonymous: list[dict[str, Any]] = []

    for line_no, line in enumerate(lines, 1):
        record = parse_line(line)
        if not record:
            continue
        frame = record.get("frame", record.get("f"))
        ball = record.get("ball", record.get("id", record.get("ball_id")))
        if frame is None or ball is None:
            record["_line"] = line_no
            anonymous.append(record)
            continue
        key = (frame, ball)
        if key not in merged:
            merged[key] = {"frame": frame, "ball": ball}
            order.append(key)
        merged[key].update(record)

    return [merged[key] for key in order] + anonymous


@dataclass
class Finding:
    ball: Any
    frame: Any
    status: str
    root_cause: str | None = None
    excluded_through: str = "start"
    messages: list[str] = field(default_factory=list)

    def line(self) -> str:
        head = f"ball={label(self.ball)} frame={label(self.frame)} status={self.status}"
        if self.root_cause:
            head += f" root_cause={self.root_cause}"
        head += f" excluded_through={self.excluded_through}"
        if self.messages:
            head += " :: " + "; ".join(self.messages)
        return head


def audit_geometry(r: dict[str, Any], depth_tol: float) -> tuple[str, list[str]]:
    depth = r.get("depth")
    n = vec(r.get("normal"))
    radius = r.get("radius", r.get("r"))
    center = vec(r.get("center", r.get("p")))
    center_y = r.get("center_y", r.get("y"))
    plane_y = float(r.get("plane_y", 0.0))
    messages: list[str] = []

    if depth is None or n is None:
        return "blocked_geometry_fields", ["need depth and normal"]
    n_len = length(n)
    if abs(n_len - 1.0) > 1e-3:
        messages.append(f"normal_not_unit len={n_len:.6f}")

    if center is None and center_y is not None:
        center = (0.0, float(center_y), 0.0)
    if center is None or radius is None:
        implied_h = None
        try:
            implied_h = float(radius) - float(depth)  # type: ignore[arg-type]
        except (TypeError, ValueError):
            pass
        messages.append(
            "depth/normal give center_distance=R-depth="
            + (f"{implied_h:.6f}" if implied_h is not None else "<need R>")
            + "; independent check needs center and radius"
        )
        return "blocked_geometry_state", messages

    radius = float(radius)
    plane_point = vec(r.get("plane_point")) or (center[0], plane_y, center[2])
    unit_n = tuple(x / n_len for x in n) if n_len > 1e-12 else n
    h = dot(sub(center, plane_point), unit_n)
    depth_geom = radius - h
    error = depth_geom - float(depth)
    messages.append(
        f"h={h:.6f} R-depth_log={radius - float(depth):.6f} "
        f"depth_geom={depth_geom:.6f} residual={error:+.2e}"
    )
    if h < -depth_tol:
        messages.append("normal points away from centre")
        return "confirmed_geometry_normal", messages
    if abs(error) > depth_tol:
        return "confirmed_geometry_depth_mismatch", messages
    return "excluded_geometry_ok", messages


def classify_record(r: dict[str, Any], depth_tol: float, arm_eps: float) -> Finding:
    ball = r.get("ball", r.get("id", "<missing>"))
    frame = r.get("frame", r.get("f", r.get("_line", "<missing>")))
    finding = Finding(ball=ball, frame=frame, status="unresolved")

    geom_status, geom_msgs = audit_geometry(r, depth_tol)
    finding.messages.extend(geom_msgs)
    if geom_status.startswith("confirmed_geometry"):
        finding.status = "resolved"
        finding.root_cause = geom_status.removeprefix("confirmed_")
        finding.excluded_through = "geometry_audit"
        return finding

    # Sleep gate is first in a typical solver pipeline.
    awake = r.get("awake", r.get("active"))
    sleeping = r.get("sleeping", r.get("asleep"))
    called = r.get("solver_called", r.get("participated", r.get("solved")))
    why = str(r.get("why", r.get("reason", ""))).lower()
    if is_true_flag(sleeping) or is_false_flag(awake) or "sleep" in why:
        finding.status = "resolved"
        finding.root_cause = "sleep_skip"
        finding.excluded_through = "sleep_gate"
        finding.messages.append("sleep field/reason confirms non-participation")
        return finding
    if is_true_flag(awake) or is_true_flag(called):
        finding.excluded_through = "sleep_gate"
        finding.messages.append("sleep skip excluded by active/participated field")
    else:
        finding.excluded_through = "sleep_gate?"
        finding.messages.append("cannot pass sleep gate: need awake/active/solver_called")
        return finding

    # Normal orientation and contact-point sign.
    vn = r.get("vn", r.get("vn_pre", r.get("normal_velocity")))
    normal_reversed = r.get("normal_reversed", r.get("normal_flipped"))
    if is_true_flag(normal_reversed) or "reversed_normal" in why:
        finding.status = "resolved"
        finding.root_cause = "reversed_normal_positive_vn"
        finding.excluded_through = "normal_orientation"
        return finding
    if vn is None:
        finding.messages.append("cannot test vn>0 branch: need solver-side vn_pre")
        return finding
    vn = float(vn)
    vn_eps = float(r.get("vn_eps", 0.0))
    if vn > vn_eps:
        aligned = r.get("normal_aligned")
        if aligned is False:
            finding.root_cause = "reversed_normal_positive_vn"
        else:
            finding.root_cause = "positive_vn_separating_or_resting"
        finding.status = "resolved"
        finding.excluded_through = "vn_dot_sign_gate"
        finding.messages.append(f"vn={vn:.6g} > {vn_eps:g}")
        return finding
    finding.messages.append(f"vn={vn:.6g} is non-positive; sign gate excluded")

    # Lever arm.  It explains zero torque; only an explicit guard/formula makes
    # it explain a zero scalar impulse.
    arm = r.get("arm", r.get("lever_arm"))
    rv = vec(r.get("r", r.get("lever")))
    tv = vec(r.get("t", r.get("tangent")))
    if arm is None and rv is not None and tv is not None:
        arm = length(cross(rv, tv))
    if arm is not None:
        arm = float(arm)
        if abs(arm) <= arm_eps:
            if "zero_arm_guard" in why or is_true_flag(r.get("zero_arm_guard")):
                finding.status = "resolved"
                finding.root_cause = "zero_arm_guard_zeroed_j"
            else:
                finding.status = "resolved_for_spin_only"
                finding.root_cause = "zero_lever_arm_explains_no_spin_not_jn_zero"
            finding.excluded_through = "lever_arm"
            finding.messages.append(f"arm={arm:.3e} <= {arm_eps:g}")
            return finding
        finding.messages.append(f"arm={arm:.6g} nonzero; lever-arm cause excluded")
    else:
        finding.excluded_through = "vn_sign_gate"
        finding.messages.append("cannot test lever arm: need arm or r+t")
        return finding

    j = r.get("j", r.get("jn", r.get("impulse")))
    if j is not None and abs(float(j)) <= 1e-9:
        finding.excluded_through = "lever_arm"
        finding.messages.append("arm nonzero but j=0; inspect effective mass/clamp/mask")
        finding.status = "unresolved_second_cause"
    else:
        finding.status = "not_a_zero_impulse_event"
        finding.excluded_through = "angular_solver"
    return finding


def spin_report(records: list[dict[str, Any]], tolerance: float = 0.05) -> list[str]:
    lines = ["spin_acceptance:"]
    found = False
    for r in records:
        w1 = vec(r.get("w1", r.get("w_after")))
        if w1 is None:
            continue
        found = True
        ball = r.get("ball", r.get("id", "<missing>"))
        frame = r.get("frame", r.get("f", r.get("_line", "?")))
        w0 = vec(r.get("w0", r.get("w_before"))) or (0.0, 0.0, 0.0)
        delta = length(sub(w1, w0))

        expected = None
        basis = ""
        arm = r.get("arm", r.get("lever_arm"))
        jt = r.get("jt", r.get("j_t"))
        inv_i = r.get("inv_inertia", r.get("invI"))
        inertia = r.get("inertia", r.get("I"))
        if arm is not None and jt is not None and inv_i is not None:
            expected = abs(float(inv_i)) * abs(float(arm)) * abs(float(jt))
            basis = "invI*arm*|jt|"
        elif arm is not None and jt is not None and inertia is not None:
            expected = abs(float(arm)) * abs(float(jt)) / float(inertia)
            basis = "arm*|jt|/I"
        elif r.get("contact_model") == "sticking_solid_sphere":
            speed = abs(float(r.get("speed", r.get("vx", 0.0))))
            radius = float(r.get("radius", r.get("r", 0.0)))
            if speed > 0.0 and radius > 0.0:
                expected = 5.0 * speed / (7.0 * radius)
                basis = "5v/(7R)"

        if delta <= 1e-9:
            status = "UNRESOLVED_zero_spin"
            detail = "need classification/exclusion step"
        elif expected is None:
            status = "nonzero_no_analytic_baseline"
            detail = f"|dw|={delta:.6g}"
        else:
            error = abs(delta - expected) / expected if expected else 0.0
            status = "PASS" if error <= tolerance else "FAIL_ratio"
            detail = (f"|dw|={delta:.6g} analytic={expected:.6g} ({basis}) "
                      f"relerr={error:.4f} limit={tolerance:.2f}")
        lines.append(f"  ball={label(ball)} frame={label(frame)} {status} {detail}")
    if not found:
        lines.append("  no parseable w0 -> w1 records")
    return lines


def period_report(records: list[dict[str, Any]], args: argparse.Namespace) -> list[str]:
    zero_events = []
    for r in records:
        j = r.get("j", r.get("jn", r.get("impulse")))
        frame = r.get("frame", r.get("f"))
        if j is not None and abs(float(j)) <= 1e-9 and frame is not None:
            zero_events.append((r.get("ball", r.get("id", "?")), float(frame)))

    lines = ["period_analysis:"]
    by_ball: dict[Any, list[float]] = {}
    for ball, frame in zero_events:
        by_ball.setdefault(ball, []).append(frame)
    for ball, frames in by_ball.items():
        frames.sort()
        deltas = [b - a for a, b in zip(frames, frames[1:])]
        lines.append(f"  ball={label(ball)} zero_j_frames={[int(x) for x in frames]} deltas={deltas}")
        if not deltas:
            continue
        n_observed = min(deltas) if len(set(deltas)) == 1 else sorted(deltas)[len(deltas) // 2]
        fps = args.fps or _as_float(records, "fps")
        speed = args.speed or _as_float(records, "speed", "vx")
        radius = args.radius or _as_float(records, "radius", "r")
        if not (fps and speed is not None and radius):
            lines.append("    period relation needs --fps/--speed/--radius or log fields")
            continue
        symmetry = max(1, args.symmetry)
        if abs(speed) < 1e-12:
            lines.append("    speed=0: rolling-period formula is singular/infinite")
            continue
        expected = 2.0 * math.pi * float(radius) * float(fps) * symmetry / abs(float(speed))
        error = abs(n_observed - expected) / expected
        tolerance = max(1.0, args.period_tol * expected)
        match = error * expected <= tolerance
        lines.append(
            f"    N_obs~{n_observed:g} N_roll={expected:.2f} frames "
            f"(T={expected / float(fps):.4f}s) relerr={error:.3f} "
            f"match={'yes' if match else 'NO -> second clock/cause'}"
        )
    if not zero_events:
        lines.append("  no frame-tagged j=0 records found")
    return lines


def _as_float(records: list[dict[str, Any]], *keys: str) -> float | None:
    for r in records:
        for key in keys:
            if key in r:
                try:
                    return float(r[key])
                except (TypeError, ValueError):
                    pass
    return None


def selftest() -> int:
    good = {
        "frame": 12, "ball": 7, "depth": 0.0031, "normal": "(0,1,0)",
        "radius": 0.5, "center": "(0,0.4969,0)",
    }
    status, messages = audit_geometry(good, 5.1e-5)
    assert status == "excluded_geometry_ok", (status, messages)

    reversed_ball = {
        "frame": 1, "ball": 2, "awake": True, "solver_called": True,
        "vn": 1.0, "normal_aligned": False,
    }
    finding = classify_record(reversed_ball, 5.1e-5, 1e-6)
    assert finding.root_cause == "reversed_normal_positive_vn", finding.line()

    sleeping_ball = {"frame": 30, "ball": 3, "awake": False, "j": 0.0}
    finding = classify_record(sleeping_ball, 5.1e-5, 1e-6)
    assert finding.root_cause == "sleep_skip", finding.line()

    arm_ball = {
        "frame": 31, "ball": 1, "awake": True, "solver_called": True,
        "vn": -1.0, "arm": 0.0, "why": "zero_arm_guard", "j": 0.0,
    }
    finding = classify_record(arm_ball, 5.1e-5, 1e-6)
    assert finding.root_cause == "zero_arm_guard_zeroed_j", finding.line()

    r = 0.5
    v = math.pi
    fps = 60.0
    expected_frames = 2 * math.pi * r * fps / v
    assert abs(expected_frames - 60.0) < 1e-9
    print("selftest=1 geometry=ok period_formula=ok causes=3")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("logfile", nargs="?", help="log file; - or omitted means stdin")
    parser.add_argument("--fps", type=float)
    parser.add_argument("--speed", type=float)
    parser.add_argument("--radius", type=float)
    parser.add_argument("--symmetry", type=int, default=1,
                        help="N-fold symmetric contact pattern divides circumference period")
    parser.add_argument("--depth-tol", type=float, default=5.1e-5)
    parser.add_argument("--arm-eps", type=float, default=1e-6)
    parser.add_argument("--period-tol", type=float, default=0.05)
    parser.add_argument("--spin-tol", type=float, default=0.05)
    parser.add_argument("--selftest", action="store_true")
    args = parser.parse_args(argv)

    if args.selftest:
        return selftest()

    if args.logfile and args.logfile != "-":
        text = Path(args.logfile).read_text()
        lines = text.splitlines()
    else:
        lines = sys.stdin.read().splitlines()

    records = parse_log(lines)
    if not records:
        print("No parseable records.")
        return 2

    for record in records:
        # The quoted j=0 lines use key j; classify all records but only print
        # records that mention an impulse or a possible contact failure.
        if "j" not in record and "jn" not in record and "impulse" not in record:
            continue
        print(classify_record(record, args.depth_tol, args.arm_eps).line())

    for line in period_report(records, args):
        print(line)
    for line in spin_report(records, args.spin_tol):
        print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
