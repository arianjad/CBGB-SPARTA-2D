#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../.." && pwd -P)
out=$(bash "$ROOT/tools/wsl/build_sparta_plain.sh" --print-plan)
grep -q '95b9abaa8bd548991cc3c3f1c58b34722f7ade74' <<< "$out"
grep -q -- '-DSPARTA_MACHINE=mpi' <<< "$out"
grep -q -- '-DSPARTA_ENABLE_TESTING=ON' <<< "$out"
grep -q 'ctest --test-dir' <<< "$out"
echo "PASS plain SPARTA build plan"
