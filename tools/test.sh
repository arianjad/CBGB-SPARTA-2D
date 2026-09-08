#!/usr/bin/env bash
# Portable checks; no solver launch, package download, or private run data.
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
cd "$repo"
python_cmd=${PYTHON:-python3}
"$python_cmd" -m compileall -q tools cases/b5-lean/gen_b5.py
"$python_cmd" tools/test_tracer_analyze.py
"$python_cmd" tools/test_tracer_common.py
"$python_cmd" tools/test_plot_trajectories.py
bash tools/wsl/test_build_sparta_plain.sh
bash cases/b5-lean/test_run_lean.sh
PYTHON="$python_cmd" bash tools/test_standalone.sh
for script in tools/*.sh tools/wsl/*.sh cases/b5-lean/*.sh; do bash -n "$script"; done
"$python_cmd" - <<'PY'
from pathlib import Path
import hashlib
import subprocess
import sys
files = sorted(Path('cases/b5-lean').glob('*.surf'))
before = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
subprocess.run([sys.executable, 'cases/b5-lean/gen_b5.py'], check=True)
assert all(hashlib.sha256(p.read_bytes()).hexdigest() == h for p, h in before.items())
print('PASS geometry regeneration matches tracked surfaces')
PY
