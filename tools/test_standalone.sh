#!/usr/bin/env bash
# A failed or unfinished helium run must not start Julia or create outputs.
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
tmp=$(mktemp -d)
trap 'rm -rf -- "$tmp"' EXIT
mkdir "$tmp/he"
touch "$tmp/he/field.grid" "$tmp/he/cell.surf"
printf '0\n' > "$tmp/he/rc.sentinel"
for status in running failed; do
  printf '{"status":"%s"}\n' "$status" > "$tmp/he/manifest.json"
  if JULIA=true bash "$repo/tools/run_2d_standalone.sh" "$tmp/he" "$tmp/output" \
      > "$tmp/log" 2>&1; then
    echo "FAIL: accepted helium status=$status" >&2; exit 1
  fi
  grep -q "helium run is not complete: status=$status" "$tmp/log"
  test ! -e "$tmp/output"
done
echo 'PASS unfinished/failed helium rejected before molecule output creation'
