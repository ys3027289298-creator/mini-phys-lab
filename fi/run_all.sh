#!/usr/bin/env bash
# One-command reproduction of the fault-injection study.
#
#   ./fi/run_all.sh
#
# Steps: build adapters -> run 9-case suite (verdicts) -> coverage/time decision.
# All outputs land in fi/results/.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"

echo "== build adapters =="
"$here/build.sh"

echo
echo "== 9-case injection suite =="
suite_rc=0
python3 "$here/tools/run_suite.py" || suite_rc=$?

echo
echo "== coverage / time decision =="
python3 "$here/tools/tradeoff.py"

echo
echo "== supplementary: shared-slot order dependence (collide) =="
"$here/build/collide" X_alias

exit "$suite_rc"
