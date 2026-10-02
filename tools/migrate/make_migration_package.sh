#!/bin/bash
# Build the ORCD -> subMIT migration package (see docs/MIGRATION-subMIT.md).
# Re-runnable: it overwrites the package, so run it again after new results land
# (e.g. the uniform-0.20mm grid-convergence tracer).
#   bash tools/migrate/make_migration_package.sh        [OUT=dir to override]
# Copies only small/derived data: frozen fields, tracer outputs, inputs, logs,
# figures. Raw field.grid dumps and restart files are left on ORCD scratch.
set -euo pipefail
R="$HOME/CBGB-SPARTA-2D"
B="$HOME/orcd/scratch/cbgb-sparta-3d"
OUT="${OUT:-$HOME/orcd/scratch/cbgb-migrate}"
mkdir -p "$OUT/claude"
cd "$R"
[ -z "$(git status --porcelain --untracked-files=no)" ] || { echo "commit tracked changes first" >&2; exit 1; }

echo "== repo bundle"
git bundle create "$OUT/CBGB-SPARTA-2D.bundle" --all
git bundle verify "$OUT/CBGB-SPARTA-2D.bundle" >/dev/null

echo "== post-processed simulation data"
dirs=""; for d in fluor-scan fluor-scan-fine fluor-scan-ports fluor-cell-3d-ne-sealed; do [ -d "$B/$d" ] && dirs="$dirs $d"; done
nice -n 19 tar -C "$B" -cf "$OUT/cbgb-sparta-3d-post.tar" \
  --exclude='field.grid' --exclude='restart.*' --exclude='*.err' --exclude='*.out.tmp' \
  --exclude='lock.d' --exclude='run.log' --exclude='moleculein.out' --exclude='moleculeout.out' \
  --exclude='.stats_jobs*' $dirs

echo "== figures + exports"
tar -C "$R" -cf "$OUT/results-figures.tar" results/figures results/exports
cp "$R/.local/julia-1.9.4-linux-x86_64.tar.gz" "$OUT/"

echo "== Claude Code files (no credentials)"
P="$HOME/.claude/projects/-home-heeric"
rm -rf "$OUT/claude/memory"; cp -r "$P/memory" "$OUT/claude/memory"
cp "$HOME/.claude/settings.json" "$OUT/claude/"
tar -C "$P" -czf "$OUT/claude/transcripts.tgz" --exclude=memory .
cp "$R/docs/MIGRATION-subMIT.md" "$OUT/README.md"

ls -la "$OUT" "$OUT/claude" | awk 'NF>=9 {printf "%10.2f MB  %s\n", $5/1e6, $9}'
echo "package ready in $OUT -- copy it with the rsync command in README.md"
