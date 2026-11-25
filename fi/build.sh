#!/usr/bin/env bash
# Build every fault-injection adapter.
#   fi/build/adapter                 plain (for correctness runs)
#   fi/build/adapter_cov             instrumented with gcov (for coverage runs)
# Adapters are wrapped so each instrumented run gets its own .gcda via
# GCOV_PREFIX, keeping per-injection-point coverage isolated.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
root="$(cd "$here/.." && pwd)"
out="$root/fi/build"
mkdir -p "$out"

cc="${CC:-gcc}"
common=(-std=c11 -Wall -Wextra -Wno-unused-function)

for m in collide solve integrate; do
  src="$root/fi/adapters/${m}_adapter.c"
  "$cc" -O2 "${common[@]}" "$src" -o "$out/$m" -lm
  "$cc" -O0 -g -fno-inline -fno-omit-frame-pointer --coverage "${common[@]}" \
        "$src" -o "$out/${m}_cov" -lm
  echo "built $m"
done
