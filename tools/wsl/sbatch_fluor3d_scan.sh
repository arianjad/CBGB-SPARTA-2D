#!/bin/bash
#SBATCH --job-name=fluor3d-scan
#SBATCH --account=mit_general
#SBATCH --partition=mit_quicktest
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=00:15:00
#SBATCH --requeue
#SBATCH --output=%x_%j.log

# Sealed fluor-cell-3d neon flow-rate scan, one flow per chain:
#   sbatch --export=ALL,SCCM=8 tools/wsl/sbatch_fluor3d_scan.sh
#
# Same case as sbatch_fluor3d_sealed.sh (the 4 SCCM reference), with:
# - MDOT, FILLN (aperture-law initial fill) and FNUM scaled linearly with
#   SCCM, which keeps ~1M simulated particles.
# - DT kept <= ~half the Ne-Ne mean collision time at the expected chamber
#   density (2.6us x 4/SCCM).
# - 0.3 s simulated per flow.
#
# The 8x8x8 chamber grid (0.41mm) resolves the mean free path only up to
# ~4-8 SCCM (lambda ~ 0.4mm x 4/SCCM). Higher flows are under-resolved:
# the stagnant-chamber state is fine, the port flow only approximate.
#
# How the chain runs:
# - Default partition is mit_quicktest. Each job runs MAXCHUNKS
#   averaging windows (sized to fit 15 min), then resubmits itself until
#   the flow reaches TOTAL. Pass CHAIN=0 to stop that; NICE=<n> lowers
#   the chain's priority; CHAIN_TIME=<hh:mm:ss> sets each job's walltime, e.g.
#   short-chunk chains on mit_normal/mit_preemptable, where short jobs backfill
#   much sooner than long ones:
#     sbatch -p mit_preemptable -t 00:30:00 --export=ALL,SCCM=10,MAXCHUNKS=4,CHAIN_TIME=00:30:00 ...
# - The same script also works on other partitions (sbatch -p ... -t ...;
#   pass MAXCHUNKS=0 for "run to the end").
# - A per-flow lock (mkdir is atomic) lets only one job at a time work on a
#   flow, whichever partition starts first; a job finding the lock held by
#   another live job re-queues itself for 30 min later, so the chain
#   survives if that other job dies (e.g. a mit_normal walltime kill).
# - EXCLUDE=<nodes> is passed to every job in the chain (slow nodes). It goes
#   through --export, which splits on commas: use a range like node[2000-2027].
# - Resumes from the newest complete checkpoint and trims field.grid back
#   to that step first, so a preempted, requeued or timed-out job never
#   leaves stray frames.
# - CHECKPOINT_STEPS = DUMPFREQ, so every restart lands on an
#   averaging-window edge.
# - NP must stay 40 for the whole flow.

set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load gcc/12.2.0 openmpi/4.1.4

: "${SCCM:?pass SCCM=<flow> via --export}"
source "$SLURM_SUBMIT_DIR/tools/wsl/fluor3d_flow_params.sh"
flow_params "$SCCM"   # DT, DUMPFREQ, CHUNKS_QT, TOTAL -- see that file
MAXCHUNKS="${MAXCHUNKS:-$CHUNKS_QT}"
read MDOT FILLN FNUM <<<"$(python3 -c "
f=$SCCM
print(f'{1.500541e-8*f:.6e} {1.35e21*f:.4e} {5e10*f:.4e}')")"

CASEDIR=cases/fluor-cell-3d
SCRIPT="$SLURM_SUBMIT_DIR/tools/wsl/sbatch_fluor3d_scan.sh"
RUNDIR="$HOME/orcd/scratch/cbgb-sparta-3d/fluor-scan/sccm$SCCM"
mkdir -p "$RUNDIR"
cp -n "$CASEDIR"/in.ne_fluor3d_sealed_run "$CASEDIR"/in.ne_fluor3d_sealed_resume \
      "$CASEDIR"/cell_fluor3d_sealed.surf "$CASEDIR"/ne.species "$CASEDIR"/ne.vss "$RUNDIR/"
cd "$RUNDIR"

resubmit() {  # $1 = extra sbatch args (e.g. a --begin delay)
  local next i
  for i in $(seq 60); do   # retry up to 30 min: quicktest caps queued jobs at 8/user
    next=$(cd "$SLURM_SUBMIT_DIR" && sbatch --parsable --partition="$SLURM_JOB_PARTITION" --nice="${NICE:-0}" \
           --time="${CHAIN_TIME:-00:15:00}" ${EXCLUDE:+--exclude=$EXCLUDE} $1 \
           --export=ALL,SCCM="$SCCM",MAXCHUNKS="$MAXCHUNKS",NICE="${NICE:-0}",CHAIN_TIME="${CHAIN_TIME:-00:15:00}",EXCLUDE="${EXCLUDE:-}" "$SCRIPT") && break
    sleep 30
  done
  echo "$(date) job $SLURM_JOB_ID: queued next ${next:-FAILED}" | tee -a "$RUNDIR/scan.log"
}

LOCK=lock.d
if ! mkdir "$LOCK" 2>/dev/null; then
  holder=$(cat "$LOCK/job" 2>/dev/null)
  if [ -n "$holder" ] && [ "$holder" != "$SLURM_JOB_ID" ] && \
     squeue -j "$holder" -h -t R 2>/dev/null | grep -q .; then
    echo "$(date) job $SLURM_JOB_ID: flow held by running job $holder" | tee -a scan.log
    [ "${CHAIN:-1}" != 0 ] && resubmit "--begin=now+1800"   # check back in 30 min
    exit 0
  fi
fi
echo "$SLURM_JOB_ID" > "$LOCK/job"
trap 'rm -rf "$RUNDIR/$LOCK"' EXIT

# newest checkpoint with all NP per-proc files + base file present and non-empty
LAST=0
for b in $(ls restart.base.* 2>/dev/null | sed 's/restart.base.//' | sort -n -r); do
  n=$(ls restart.[0-9]*."$b" 2>/dev/null | wc -l)
  if [ "$n" -eq "$SLURM_NTASKS" ] && [ -z "$(find . -maxdepth 1 -name "restart.*.$b" -size 0)" ]; then
    LAST=$b; break
  fi
done
# keep only the two newest checkpoints
for b in $(ls restart.base.* 2>/dev/null | sed 's/restart.base.//' | sort -n -r | tail -n +3); do
  [ "$b" -lt "$LAST" ] && rm -f restart.*."$b"
done
REMAIN=$(( (TOTAL - LAST) / DUMPFREQ ))
if [ "$REMAIN" -le 0 ]; then echo "already complete at step $LAST" | tee -a scan.log; exit 0; fi
N=$REMAIN; [ "$MAXCHUNKS" -gt 0 ] && [ "$MAXCHUNKS" -lt "$N" ] && N=$MAXCHUNKS
echo "SCCM=$SCCM MDOT=$MDOT FILLN=$FILLN FNUM=$FNUM DT=$DT DUMPFREQ=$DUMPFREQ TOTAL=$TOTAL" >> scan.log
COMMON="-var MDOT $MDOT -var FILLN $FILLN -var FNUM $FNUM -var DT $DT \
        -var DUMPFREQ $DUMPFREQ -var CHECKPOINT_STEPS $DUMPFREQ -var NCHUNKS $N"

if [ "$LAST" -eq 0 ]; then
  echo "$(date) job $SLURM_JOB_ID ($SLURM_JOB_PARTITION): fresh start, $N of $REMAIN chunks" | tee -a scan.log
  mpirun -np "$SLURM_NTASKS" "$HOME/opt/sparta-27Aug2026/bin/spa_mpi" \
    -in in.ne_fluor3d_sealed_run $COMMON &> run.log
else
  # drop any frames dumped after the checkpoint (including a truncated one)
  if [ -f field.grid ]; then
    cut=$(grep -b -A1 "^ITEM: TIMESTEP" field.grid | awk -F'[-:]' -v s="$LAST" \
          'prev != "" && $2+0 > s {print prev; exit} /ITEM: TIMESTEP/ {prev=$1; next} {prev=""}')
    [ -n "$cut" ] && truncate -s "$cut" field.grid
  fi
  echo "$(date) job $SLURM_JOB_ID ($SLURM_JOB_PARTITION): resume from step $LAST, $N of $REMAIN chunks" | tee -a scan.log
  mpirun -np "$SLURM_NTASKS" "$HOME/opt/sparta-27Aug2026/bin/spa_mpi" \
    -in in.ne_fluor3d_sealed_resume -var RESTARTFILE "restart.%.$LAST" $COMMON &>> run.log
fi
rc=$?
echo "$(date) job $SLURM_JOB_ID: solver exit code $rc" | tee -a scan.log
tail -2 run.log
if [ "$rc" -eq 0 ] && [ "${CHAIN:-1}" != 0 ] && [ $(( LAST + N * DUMPFREQ )) -lt "$TOTAL" ]; then
  rm -rf "$RUNDIR/$LOCK"; trap - EXIT
  resubmit ""
fi
exit "$rc"
