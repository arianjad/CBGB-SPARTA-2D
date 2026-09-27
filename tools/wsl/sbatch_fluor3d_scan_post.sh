#!/bin/bash
#SBATCH --job-name=fluor3d-scan-post
#SBATCH --account=mit_general
#SBATCH --partition=mit_quicktest
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=00:15:00
#SBATCH --output=%x_%j.log

# Post-processing for one flow of the fluor-cell-3d scan (sbatch_fluor3d_scan.sh):
#   sbatch --export=ALL,SCCM=8[,NPART=1000] tools/wsl/sbatch_fluor3d_scan_post.sh
# Writes the following into RUNDIR/post/:
#   1. asymmetry.txt: half-space asymmetry and settling check over the last
#      75 ms (tools/fluor3d_asymmetry.py).
#   2. frozen.grid: the last 20 ms time-averaged into one frame
#      (tools/field3d_freeze_avg.py).
#   3. tracer_<cond>.out and tracer_summary.txt: the BaF+ tracer at the
#      three handoff cross-section conditions (energy-dependent default, and
#      constant 2.0e-18 and 1.5e-19), classified by tools/tracer_classify3d.py.
# 8 cores so it fits beside a running 40-core scan job under quicktest's
# 48-core/user cap. Tracer cost grows with chamber density (more collisions
# to diffuse out), so pass a smaller NPART, or a longer partition/time, at
# high flow, or NPART=0 to only prepare the field.

set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
: "${SCCM:?pass SCCM=<flow> via --export}"
NPART="${NPART:-1000}"
# DT and averaging-window length (steps): shared with sbatch_fluor3d_scan.sh
source tools/wsl/fluor3d_flow_params.sh
flow_params "$SCCM"; WIN=$DUMPFREQ
# freeze: last 20 ms; asymmetry: last 75 ms (read from the dump's tail only)
read SPAN NFREEZE NASYM <<<"$(python3 -c "
import math; dt=$DT; w=$WIN
print(int(round(0.02/dt)), max(1, math.ceil(0.02/dt/w)), max(1, math.ceil(0.075/dt/w)))")"
NASYM="${ASYM_FRAMES:-$NASYM}"   # override for huge dumps (127 SCCM: 32 frames x 340MB)
RUNDIR="$HOME/orcd/scratch/cbgb-sparta-3d/fluor-scan/sccm$SCCM"
[ "$SCCM" = 4 ] && RUNDIR="$HOME/orcd/scratch/cbgb-sparta-3d/fluor-cell-3d-ne-sealed"
OUT="$RUNDIR/post"; mkdir -p "$OUT"

module purge
module load deprecated-modules gcc/12.2.0-x86_64 python/3.10.8-x86_64
source .venv/bin/activate
python tools/fluor3d_asymmetry.py "$RUNDIR/field.grid" "$NASYM" > "$OUT/asymmetry.txt"
python tools/field3d_freeze_avg.py "$RUNDIR/field.grid" --out "$OUT/frozen.grid" --span-steps "$SPAN" --max-frames "$NFREEZE"
cp -n cases/fluor-cell-3d/cell_fluor3d_sealed.surf "$OUT/"

# NPART=0: field prep only -- run the tracer afterwards in pieces with
# sbatch_fluor3d_tracer_chunk.sh (needed at high flow, where 3 x 1000 ions
# don't fit in 15 min)
[ "$NPART" = 0 ] && { echo "prep only (NPART=0)"; exit 0; }
JULIA="$SLURM_SUBMIT_DIR/.local/julia-1.9.4/bin/julia"
cd "$OUT"
for C in "edep:" "const2e-18:--sigma-energy-dependent 0 --sigma 2.0e-18" \
         "const1.5e-19:--sigma-energy-dependent 0 --sigma 1.5e-19"; do
  N=${C%%:*}; A=${C#*:}
  "$JULIA" --project="$SLURM_SUBMIT_DIR/tracer" --threads="$SLURM_CPUS_PER_TASK" \
    "$SLURM_SUBMIT_DIR/tracer/ParticleTracing3D.jl" cell_fluor3d_sealed.surf frozen.grid \
    -n "$NPART" -x 0.03048 -y 0.01905 -z 0.01905 --vx 24843.71 --vy 0 --vz 0 \
    -m 20.1797 -M 156.325 $A --saveall 1 --seed 42 > "tracer_$N.out" 2> "tracer_$N.err"
  echo "$N julia exit $?"
done
cd "$SLURM_SUBMIT_DIR"
python tools/tracer_classify3d.py "$OUT"/tracer_edep.out "$OUT"/tracer_const2e-18.out \
  "$OUT"/tracer_const1.5e-19.out | tee "$OUT/tracer_summary.txt"
