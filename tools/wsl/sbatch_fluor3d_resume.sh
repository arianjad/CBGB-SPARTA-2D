#!/bin/bash
#SBATCH --job-name=fluor3d-resume
#SBATCH --partition=mit_quicktest
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=00:15:00
#SBATCH --output=%x_%j.log

set -uo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load gcc/12.2.0 openmpi/4.1.4

CASEDIR=cases/fluor-cell-3d
RUNDIR="$HOME/orcd/scratch/cbgb-sparta-3d/fluor-cell-3d-ne"
cp -n "$CASEDIR"/in.ne_fluor3d_resume "$RUNDIR/"
cd "$RUNDIR"

# mit_quicktest: mit_normal had ~1100-1400 jobs queued, and 6 chunks x 20000
# steps took 9m32s, inside the 15-min cap. If it hits the cap, just resume
# again from the last restart.*.<step> file.
#
# Keep NP=40 (same as the run that wrote RESTARTFILE). SPARTA can redistribute
# to a different count on read, but a mismatched-NP test here failed with
# "ghost cells do not exist".
RESTARTFILE="${RESTARTFILE:-restart.%.120000}"
NCHUNKS="${NCHUNKS:-6}"
CHECKPOINT_STEPS="${CHECKPOINT_STEPS:-20000}"

mpirun -np "$SLURM_NTASKS" "$HOME/opt/sparta-27Aug2026/bin/spa_mpi" \
  -in in.ne_fluor3d_resume \
  -var RESTARTFILE "$RESTARTFILE" \
  -var MDOT 6.002164e-08 -var FILLN 5.4e21 -var FNUM 1e11 -var DT 1e-7 \
  -var DUMPFREQ 20000 -var CHECKPOINT_STEPS "$CHECKPOINT_STEPS" -var NCHUNKS "$NCHUNKS" \
  &> resume.log
rc=$?
echo "solver exit code: $rc"
tail -30 resume.log
exit "$rc"
