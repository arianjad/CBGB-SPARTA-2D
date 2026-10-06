#!/bin/bash
#SBATCH --job-name=cyl3d-tracer
#SBATCH --partition=submit
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/%x_%j.log

# BaF+ tracer on one frozen flow of the plain-cylinder scan, at a list of
# beam energies (one Julia run each, sequentially):
#   sbatch --export=ALL,SCCM=20,ENERGIES="200 400 600",NPART=2000,SEED=1 tools/wsl/sbatch_cyl3d_tracer.sh
# Beam along +x on the axis, starting at the inlet aperture's outer face
# (x = -1.27 mm), monoenergetic (-T 0) like the fluor-cell-3d beam runs.
# COND: table-central (default), table-low, table-high or edep.
# Writes into RUNDIR/post/, per energy (tag = tracer_<COND>-e<E>_s<SEED>):
#   <tag>.out        final state of every ion (idx x y z xn yn zn vx vy vz collides time)
#   <tag>_ions.csv   per-ion milestones: where/when/after how many collisions the
#                    ion first drops below 1 eV and 0.1 eV, max depth and radius,
#                    path length, final KE (--ionlog)
#   <tag>_xr.csv     (x, r) histograms of energy deposited at collisions (-> dE/dx,
#                    the Bragg curve), ion dwell time and collisions (--hist)
#   <tag>_xe.csv     collisions binned by x and log10 E (slowing-down band)
# An energy whose .out already exists is skipped, so a timed-out job can just
# be resubmitted.
# Seeds must differ between chunks of the same (flow, COND, energy).
# GIF=1 (default): per energy, also trace 40 ions with every flight leg recorded
# (--trajprint 40, seed 300) and render post/gif/sccm<F>_<COND>_e<E>.gif (+ _paths.png,
# .npz) with tools/render_fluor3d_trajectories.py: 160 frames, log time to 1 ms,
# then linear. The raw leg file (up to ~1 GB) is deleted after rendering.

set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
: "${SCCM:?}"; : "${ENERGIES:?}"; : "${SEED:?}"
NPART="${NPART:-2000}"; COND="${COND:-table-central}"
case "$COND" in
  table-central|table-low|table-high) A="--sigma-model table --sigma-table-variant ${COND#table-}" ;;
  edep) A="--sigma-model edep" ;;
  *) echo "unknown COND $COND"; exit 2 ;;
esac
source "$SLURM_SUBMIT_DIR/tools/wsl/cyl3d_flow_params.sh"; cyl_flow_params "$SCCM"   # SURF, RUNBASE_DEF, GAS, MGAS_U for CONFIG
# helium (CONFIG=ports-he): BaF+-He table (tools/cross_section/ion_he.py), same variant column
[ "$GAS" = He ] && A="$A --sigma-table $SLURM_SUBMIT_DIR/tracer/data/bafplus_he_sigma_tracer.tsv"
GEOM=$(case "$CONFIG" in ports|ports-he) echo cylports ;; *) echo cyl ;; esac)
RUNDIR="$HOME/orcd/scratch/cbgb-sparta-3d/${RUNBASE:-$RUNBASE_DEF}/sccm$SCCM"
cd "$RUNDIR/post" || exit 2
[ -f frozen.grid ] || { echo "no frozen.grid -- run sbatch_cyl3d_post.sh first"; exit 2; }
JULIA="$SLURM_SUBMIT_DIR/.local/julia-1.9.4/bin/julia"
# Per-CPU-type Julia compile cache. The shared ~/.julia/compiled lives on NFS; a job
# on a new CPU type rewrites its .so files in place, and running jobs on other nodes
# that have them memory-mapped then die with SIGBUS ("Bus error", rc 135) hours later
# (2026-10-04: two 30/40 SCCM tracer jobs). New caches now go to a depot per CPU model.
CPUKEY=$(grep -m1 "model name" /proc/cpuinfo | md5sum | cut -c1-8)
export JULIA_DEPOT_PATH="$HOME/.julia_bycpu/$CPUKEY:$HOME/.julia"
# Free-flight cap (tracer --max-flight, null steps): without it each flight uses the gas
# density at its start point, which pushes ions from the low-density inlet bore 2-4 mm too
# deep (2026-10-05, tools/sanity/). The cyl-ports campaign was started without it and keeps
# it off (MAXFLIGHT=0) so its 200 points stay consistent, and so does the helium comparison
# run cyl-ports-he (user's choice, 2026-10-05); everything else gets 0.2 mm.
RB="${RUNBASE:-$RUNBASE_DEF}"
MAXFLIGHT="${MAXFLIGHT:-$(case "$RB" in cyl-ports|cyl-ports-he) echo 0 ;; *) echo 2e-4 ;; esac)}"
echo "RUNBASE $RB, MAXFLIGHT $MAXFLIGHT m"
fail=0
for E in $ENERGIES; do
  tag="tracer_${COND}-e${E}_s$SEED"
  [ -f "$tag.out" ] && { echo "$tag.out exists, skipping"; continue; }
  VX=$(awk -v E="$E" 'BEGIN{printf "%.2f", sqrt(2*E*1.602177e-19/(156.325*1.66054e-27))}')
  for attempt in 1 2 3; do   # a fast failure is usually a precompile collision
    t0=$(date +%s)
    "$JULIA" --project="$SLURM_SUBMIT_DIR/tracer" --threads="$SLURM_CPUS_PER_TASK" \
      "$SLURM_SUBMIT_DIR/tracer/ParticleTracing3D.jl" "$SURF" frozen.grid \
      -n "$NPART" -x -0.00127 -y 0 -z 0 --vx "$VX" --vy 0 --vz 0 -T 0 \
      -m $MGAS_U -M 156.325 $A --saveall 1 --seed "$SEED" --hist "$tag" --ionlog "${tag}_ions.csv" --max-flight "$MAXFLIGHT" \
      > "$tag.out.tmp" 2> "$tag.err"
    rc=$?
    [ "$rc" -eq 0 ] || [ $(( $(date +%s) - t0 )) -gt 300 ] && break
    echo "E=$E attempt $attempt failed fast (rc $rc), retrying"; sleep $((30 * attempt))
  done
  if [ "$rc" -eq 0 ]; then mv "$tag.out.tmp" "$tag.out"; else fail=1; fi
  echo "$(date) SCCM=$SCCM E=$E eV NPART=$NPART rc=$rc $(( $(date +%s) - t0 )) s"
done

[ "${GIF:-1}" = 1 ] || exit $fail
mkdir -p gif
source "$SLURM_SUBMIT_DIR/tools/wsl/site_env.sh"; site_python; source "$SLURM_SUBMIT_DIR/.venv/bin/activate"
export PYVISTA_OFF_SCREEN=true
for E in $ENERGIES; do
  g="gif/sccm${SCCM}_${COND}_e$E"
  [ -f "$g.gif" ] && continue
  VX=$(awk -v E="$E" 'BEGIN{printf "%.2f", sqrt(2*E*1.602177e-19/(156.325*1.66054e-27))}')
  t0=$(date +%s)
  "$JULIA" --project="$SLURM_SUBMIT_DIR/tracer" --threads="$SLURM_CPUS_PER_TASK" \
    "$SLURM_SUBMIT_DIR/tracer/ParticleTracing3D.jl" "$SURF" frozen.grid \
    -n 40 -x -0.00127 -y 0 -z 0 --vx "$VX" --vy 0 --vz 0 -T 0 \
    -m $MGAS_U -M 156.325 $A --saveall 1 --seed 300 --trajprint 40 --max-flight "$MAXFLIGHT" > "$g.traj" 2> "$g.err" \
  && python "$SLURM_SUBMIT_DIR/tools/render_fluor3d_trajectories.py" "$g.traj" --surf "$SURF" \
       --out "$g" --geom "$GEOM" --sccm "$SCCM" --energy "$E eV" --linear-after 1e-3 > "$g.render.log" 2>&1 \
  || fail=1
  rm -f "$g.traj"
  echo "$(date) SCCM=$SCCM E=$E eV gif rc=$fail $(( $(date +%s) - t0 )) s"
done
exit $fail
