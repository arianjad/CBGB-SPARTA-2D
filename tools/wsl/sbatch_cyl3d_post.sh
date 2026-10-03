#!/bin/bash
#SBATCH --job-name=cyl3d-post
#SBATCH --partition=submit
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.log

# Freeze one flow of the plain-cylinder neon scan (sbatch_cyl3d_neon.sh) into
# RUNDIR/post/frozen.grid, the last 20 ms averaged (tools/field3d_freeze_stream.py,
# one frame in memory), print the chamber summary into post/chamber.txt, and
# check the grid rule (Kn_cell >= KNMIN, DT < tau in every cell) and the
# chamber density vs the aperture law into post/rule_check.txt
# (tools/cyl3d_check_rule.py; the job fails if the rule is violated):
#   sbatch --export=ALL,SCCM=20[,CONFIG=ports|ap5][,RUNBASE=...] tools/wsl/sbatch_cyl3d_post.sh
# SPAN_S overrides the 20 ms (e.g. a short smoke run with one frame).

set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
: "${SCCM:?pass SCCM=<flow> via --export}"
source tools/wsl/cyl3d_flow_params.sh
cyl_flow_params "$SCCM"
read SPAN NFRAMES <<<"$(python3 -c "
import math; s=${SPAN_S:-0.02}; print(int(round(s/$DT)), max(1, math.ceil(s/$DT/$WINDOW)))")"
RUNDIR="$HOME/orcd/scratch/cbgb-sparta-3d/${RUNBASE:-$RUNBASE_DEF}/sccm$SCCM"
OUT="$RUNDIR/post"; mkdir -p "$OUT"
source tools/wsl/site_env.sh; site_python
source .venv/bin/activate
python tools/field3d_freeze_stream.py "$RUNDIR/field.grid" --out "$OUT/frozen.grid" \
  --span-steps "$SPAN" --max-frames "$NFRAMES" --chamber-cyl 0 45.47 15.24 > "$OUT/freeze.log"
rc=$?
grep "=== chamber" "$OUT/freeze.log" > "$OUT/chamber.txt"
cp -n "$RUNDIR/$SURF" "$OUT/"
cat "$OUT/freeze.log"
[ "$rc" -eq 0 ] || exit $rc
python tools/cyl3d_check_rule.py "$OUT/frozen.grid" --vss "$RUNDIR/ne.vss" --dt "$DT" --sccm "$SCCM" --cellfrac "$RUNDIR/cellfrac.grid" --knmin "$KN" --aeff-mm2 "$AEFF_MM2" | tee "$OUT/rule_check.txt"
exit ${PIPESTATUS[0]}
