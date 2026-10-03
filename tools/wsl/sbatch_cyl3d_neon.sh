#!/bin/bash
#SBATCH --job-name=cyl3d-neon
#SBATCH --partition=submit
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --time=2-00:00:00
#SBATCH --requeue
#SBATCH --output=logs/%x_%j.log

# Neon DSMC for one flow of the plain-cylinder cell (cases/cyl-3d), run to the
# end in one job (subMIT allows 6 days). Job size per flow comes from
# cyl3d_flow_params.sh (submit_cyl3d_campaign.sh passes it):
#   sbatch -n <NP> --mem=<MEM_GB>G -t <TIME> --export=ALL,SCCM=20 tools/wsl/sbatch_cyl3d_neon.sh
# Numerics per flow (strict cell < mean free path, dt < collision time rule,
# grid adapted at the start): tools/wsl/cyl3d_flow_params.sh. Options (env):
#   CONFIG   port geometry: ports (default) or ap5 (see cyl3d_flow_params.sh)
#   RUNBASE  data dir under ~/orcd/scratch/cbgb-sparta-3d (default cyl-ports / cyl-scan)
#   SIM_S    cap on simulated seconds; runs stop at steady state (default cap 0.5 s ports, 0.1 ap5)
#   STEPS    adapt, then run only this many production steps (smoke tests)
#   NDUMP    windows at the end that write field.grid frames (default 3)
#   KNMIN    cell <= lambda / KNMIN (default 0 = fixed lattice for ports, 0.5 for ap5)
#   DT_SET, FNUM_SET, WINDOW_SET, PPC, NADAPT  overrides
# Restartable: a requeued or resubmitted job resumes from the newest complete
# checkpoint (one per window) and trims field.grid back to it. NP must stay
# the same for the life of a flow (restart files are per processor). A
# per-flow lock keeps two jobs off one flow.

set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
source tools/wsl/site_env.sh; site_mpi

: "${SCCM:?pass SCCM=<flow> via --export}"
source tools/wsl/cyl3d_flow_params.sh
cyl_flow_params "$SCCM"   # MDOT FILLN FNUM DT REFINE ADAPT_STEPS NADAPT WINDOW TOTAL
RUNBASE="${RUNBASE:-$RUNBASE_DEF}"
CASEDIR=cases/cyl-3d
RUNDIR="$HOME/orcd/scratch/cbgb-sparta-3d/$RUNBASE/sccm$SCCM"
mkdir -p "$RUNDIR"
cp -n "$CASEDIR"/in.ne_cyl3d_run "$CASEDIR"/in.ne_cyl3d_resume "$CASEDIR"/$SURF \
      "$CASEDIR"/ne.species "$CASEDIR"/ne.vss "$RUNDIR/"
cd "$RUNDIR"
SPA="$HOME/opt/sparta-27Aug2026/bin/spa_mpi"

if [ -n "${STEPS:-}" ]; then   # smoke test: windows of STEPS steps; loose TOL forces the steady path after 2, then 1 averaging window
  echo "SMOKE SCCM=$SCCM MDOT=$MDOT FILLN=$FILLN FNUM=$FNUM DT=$DT REFINE=$REFINE ADAPT_STEPS=$ADAPT_STEPS STEPS=$STEPS NP=$SLURM_NTASKS" | tee -a scan.log
  mpirun -np "$SLURM_NTASKS" "$SPA" -in in.ne_cyl3d_run -var MDOT $MDOT -var FILLN $FILLN -var FNUM $FNUM \
    -var DT $DT -var REFINE $REFINE -var ADAPT_STEPS $ADAPT_STEPS -var NADAPT $NADAPT -var KNMIN $KN -var SURF $SURF -var XOUT $XOUT \
    -var DUMPFREQ $STEPS -var CHECKPOINT_STEPS $STEPS -var MAXWIN 3 -var MINWIN 2 -var TOL 2 -var NDUMPWIN 1 &> run.log
  rc=$?; grep -h "^ADAPT\|^STRICT RULE" log.sparta | tee -a scan.log
  grep -q "^STRICT RULE VIOLATED" log.sparta && rc=3
  echo "$(date) smoke exit $rc" | tee -a scan.log; exit $rc
fi

LOCK=lock.d
if ! mkdir "$LOCK" 2>/dev/null; then
  holder=$(cat "$LOCK/job" 2>/dev/null)
  if [ -n "$holder" ] && [ "$holder" != "$SLURM_JOB_ID" ] && squeue -j "$holder" -h -t R 2>/dev/null | grep -q .; then
    echo "$(date) job $SLURM_JOB_ID: flow held by running job $holder" | tee -a scan.log; exit 0
  fi
fi
echo "$SLURM_JOB_ID" > "$LOCK/job"
trap 'rm -rf "$RUNDIR/$LOCK"' EXIT

# newest checkpoint with all NP per-proc files present and non-empty
LAST=0
for b in $(ls restart.base.* 2>/dev/null | sed 's/restart.base.//' | sort -n -r); do
  n=$(ls restart.[0-9]*."$b" 2>/dev/null | wc -l)
  if [ "$n" -eq "$SLURM_NTASKS" ] && [ -z "$(find . -maxdepth 1 -name "restart.*.$b" -size 0)" ]; then
    LAST=$b; break
  fi
done
for b in $(ls restart.base.* 2>/dev/null | sed 's/restart.base.//' | sort -n -r | tail -n +3); do
  [ "$b" -lt "$LAST" ] && rm -f restart.*."$b"
done
# Production runs to steady state (see in.ne_cyl3d_run), capped at TOTAL steps
# (SIM_S); steady.txt records STEADY / PRODUCTION DONE across restarts.
NDUMP="${NDUMP:-3}"; MAXWIN=$(( TOTAL / WINDOW )); STEADYDONE=0; NDUMPWIN=$NDUMP; MINWIN=5
if [ "$LAST" -gt 0 ] && grep -q "PRODUCTION DONE" steady.txt 2>/dev/null; then
  echo "already complete at step $LAST" | tee -a scan.log; exit 0
fi
if [ "$LAST" -gt 0 ]; then
  S=$(grep -o "STEADY at step [0-9]*" steady.txt 2>/dev/null | tail -1 | awk '{print $4}')
  if [ -n "$S" ]; then
    STEADYDONE=1; NDUMPWIN=$(( NDUMP - (LAST - S) / WINDOW ))
    [ "$NDUMPWIN" -le 0 ] && { echo "averaging windows already done at step $LAST" | tee -a scan.log; echo "PRODUCTION DONE" >> steady.txt; exit 0; }
  else
    MAXWIN=$(( MAXWIN - LAST / WINDOW )); [ "$MAXWIN" -lt 1 ] && MAXWIN=1
    MINWIN=$(( 5 - LAST / WINDOW )); [ "$MINWIN" -lt 3 ] && MINWIN=3
  fi
fi
echo "CONFIG=$CONFIG SCCM=$SCCM MDOT=$MDOT FILLN=$FILLN FNUM=$FNUM DT=$DT REFINE=$REFINE ($CELL_MM mm) ADAPT_STEPS=$ADAPT_STEPS NADAPT=$NADAPT WINDOW=$WINDOW MAXSTEPS=$TOTAL NP=$SLURM_NTASKS" >> scan.log
COMMON="-var MDOT $MDOT -var FILLN $FILLN -var FNUM $FNUM -var DT $DT \
        -var DUMPFREQ $WINDOW -var CHECKPOINT_STEPS $WINDOW -var MAXWIN $MAXWIN -var MINWIN $MINWIN \
        -var NDUMPWIN $NDUMPWIN -var STEADYDONE $STEADYDONE"
if [ "$LAST" -eq 0 ]; then
  # (a run killed during the short adaptation phase just starts over)
  echo "$(date) job $SLURM_JOB_ID: fresh start (adapt, then up to $MAXWIN windows to steady state + $NDUMPWIN averaging)" | tee -a scan.log
  rm -f field.grid steady.txt
  mpirun -np "$SLURM_NTASKS" "$SPA" -in in.ne_cyl3d_run $COMMON \
    -var REFINE $REFINE -var ADAPT_STEPS $ADAPT_STEPS -var NADAPT $NADAPT -var KNMIN $KN -var SURF $SURF -var XOUT $XOUT &> run.log
else
  if [ -f field.grid ]; then   # drop frames dumped after the checkpoint
    cut=$(grep -b -A1 "^ITEM: TIMESTEP" field.grid | awk -F'[-:]' -v s="$LAST" \
          'prev != "" && $2+0 > s {print prev; exit} /ITEM: TIMESTEP/ {prev=$1; next} {prev=""}')
    [ -n "$cut" ] && truncate -s "$cut" field.grid
  fi
  echo "$(date) job $SLURM_JOB_ID: resume from step $LAST (steady reached: $STEADYDONE; max $MAXWIN settle windows, $NDUMPWIN averaging)" | tee -a scan.log
  mpirun -np "$SLURM_NTASKS" "$SPA" -in in.ne_cyl3d_resume -var RESTARTFILE "restart.%.$LAST" $COMMON &>> run.log
fi
rc=$?
grep -h "^ADAPT\|^STRICT RULE\|STEADY\|PRODUCTION DONE" log.sparta | tee -a scan.log
grep -q "^STRICT RULE VIOLATED" log.sparta && rc=3   # SPARTA's quit exits 0
[ "$rc" -eq 0 ] && ! grep -q "PRODUCTION DONE" steady.txt 2>/dev/null && { echo "solver ended before PRODUCTION DONE" | tee -a scan.log; rc=4; }
echo "$(date) job $SLURM_JOB_ID: solver exit code $rc" | tee -a scan.log
exit "$rc"
