#!/usr/bin/env bash
# Exercise live-field selection through the real wrapper; stop at Julia.
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
python_cmd=${PYTHON:-python3}
tmp=$(mktemp -d)
trap 'rm -rf -- "$tmp"' EXIT
mkdir "$tmp/he"
cp "$repo/cases/b5-lean/cell_b5.surf" "$tmp/he/"
"$python_cmd" - "$tmp/he/field.grid" <<'PY'
from pathlib import Path
import sys
rows = []
for i in range(10):
    for j in range(10):
        x, y = .02 + i*.001, j*.001
        rows.append(f'{10*i+j+1} {x+.0005} {y+.0005} {x} {y} {x+.001} {y+.001} 1e20 {10+i} {j} {4+i*.01}\n')
def frame(step):
    return (f'ITEM: TIMESTEP\n{step}\nITEM: NUMBER OF CELLS\n100\n'
            'ITEM: BOX BOUNDS oo ao pp\n0 .08\n0 .025\n-.001 .001\n'
            'ITEM: CELLS id xc yc xlo ylo xhi yhi nrho u v temp\n' + ''.join(rows))
Path(sys.argv[1]).write_text(frame(20000) + frame(40000) + frame(60000)[:180])
PY
cat > "$tmp/julia" <<'SH'
#!/usr/bin/env bash
if [ "$1" = --version ]; then echo 'Julia test boundary'; exit 0; fi
exit 77
SH
chmod +x "$tmp/julia"
for status in running failed; do
  printf '{"status":"%s","command":["spa_mpi","-in","in.he","-var","DT","1e-7"],"runtime_controls":{"mpi_ranks":1}}\n' "$status" > "$tmp/he/manifest.json"
  rc=0
  UNTIL_MS=5 JULIA="$tmp/julia" PYTHON="$python_cmd" \
    bash "$repo/tools/run_2d_standalone.sh" "$tmp/he" "$tmp/$status" > "$tmp/log" 2>&1 || rc=$?
  if [ "$rc" != 77 ]; then
    cat "$tmp/log" "$tmp/$status/convert.stderr" "$tmp/$status/helium-plot.stderr"
    echo "FAIL: did not reach Julia for $status" >&2; exit 1
  fi
  "$python_cmd" - "$tmp/$status" "$status" <<'PY'
from pathlib import Path
import hashlib, json, sys
p = Path(sys.argv[1])
m = json.loads((p / 'manifest.json').read_text())
assert m['current_stage'] == 'trace_crossings'
assert m['helium_manifest_status'] == sys.argv[2]
assert m['helium_snapshot']['source_timestep'] == 40000
raw = p / 'field/field.grid'
assert raw.read_text().count('ITEM: TIMESTEP') == 1
assert m['input_sha256'][str(raw)] == hashlib.sha256(raw.read_bytes()).hexdigest()
assert list((p / 'figures').glob('b5-2d-fields-*.png'))
PY
done
echo 'PASS running/failed parent: selected complete 4 ms frame, plotted frozen copy, reached Julia'
