#!/bin/bash
#SBATCH --job-name=fluor3d-run
#SBATCH --partition=mit_normal
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --output=%x_%j.log

set -uo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load gcc/12.2.0 openmpi/4.1.4

CASEDIR=cases/fluor-cell-3d
# big dump goes on scratch, not $HOME (200GB, shared)
RUNDIR="$HOME/orcd/scratch/cbgb-sparta-3d/fluor-cell-3d-ne"
PLOT_MODULE_PY="deprecated-modules gcc/12.2.0-x86_64 python/3.10.8-x86_64"

mkdir -p "$RUNDIR" results/figures
cp -n "$CASEDIR"/in.ne_fluor3d_run "$CASEDIR"/in.ne_fluor3d_resume \
      "$CASEDIR"/cell_fluor3d.surf "$CASEDIR"/ne.species "$CASEDIR"/ne.vss "$RUNDIR/"
cd "$RUNDIR"

# Restart file every CHECKPOINT_STEPS (as in cases/takahashi-3d); resume a killed
# job from the latest restart.<step> with in.ne_fluor3d_resume
# (docs/hpc-fluor-cell-3d-handoff.md). Chunks match the field.grid cadence
# (20000 steps/frame): 6 x 20000 = 120000 steps (12 ms), taken from the He case
# and not yet checked for this geometry. If not settled, resume with more chunks.
NCHUNKS=6
CHECKPOINT_STEPS=20000

# background progress plot; reads only the tail of field.grid, safe while running
(
  cd "$SLURM_SUBMIT_DIR"
  module purge
  module load $PLOT_MODULE_PY
  source .venv/bin/activate
  cd "$RUNDIR"
  while true; do
    sleep 300
    if [ ! -f field.grid ]; then continue; fi
    python "$SLURM_SUBMIT_DIR/tools/plot_field_3d.py" . \
      --out "$SLURM_SUBMIT_DIR/results/figures/fluor3d_ne_progress.png" \
      --frac 0.15 --nr 80 2>>monitor.log || true
  done
) &
monitor_pid=$!

# FNUM=1e11. The 2.5e17 from the 2D b5-lean case gives zero particles here (fill
# and mflow injection both round to 0). Quicktest scan of 1e12/5e11/1e11 (10000
# steps each) all had enough per-cell statistics and were cheap, so 1e11
# (~667 particles/occupied cell, ~30 min extrapolated) for better statistics
# near the small, under-resolved ports.
mpirun -np "$SLURM_NTASKS" "$HOME/opt/sparta-27Aug2026/bin/spa_mpi" \
  -in in.ne_fluor3d_run \
  -var MDOT 6.002164e-08 -var FILLN 5.4e21 -var FNUM 1e11 -var DT 1e-7 \
  -var DUMPFREQ 20000 -var CHECKPOINT_STEPS $CHECKPOINT_STEPS -var NCHUNKS $NCHUNKS \
  &> run.log
rc=$?
kill "$monitor_pid" 2>/dev/null
echo "solver exit code: $rc"
tail -30 run.log
exit "$rc"
