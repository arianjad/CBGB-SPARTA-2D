#!/bin/bash
#SBATCH --job-name=fluor3d-tracer
#SBATCH --account=mit_general
#SBATCH --partition=mit_quicktest
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=00:15:00
#SBATCH --output=%x_%j.log

# Extra BaF+ tracer statistics for one flow of the fluor-cell-3d scan, on
# the frozen field written by sbatch_fluor3d_scan_post.sh:
#   sbatch --export=ALL,SCCM=8,SEED=43,NPART=2000[,COND=edep] tools/wsl/sbatch_fluor3d_tracer_chunk.sh
# Writes RUNDIR/post/tracer_<COND>_s<SEED>.out. Chunks with different seeds
# are independent samples; tools/plot_fluor3d_scan.py pools every
# tracer_<COND>*.out it finds.
# COND is edep (default), table-central, table-low, table-high, const2e-18,
# const1.5e-19, or thermal-const1.93e-18 (thermal-injection test).

set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
: "${SCCM:?}"; : "${SEED:?}"; NPART="${NPART:-2000}"; COND="${COND:-edep}"
# injection: 500 eV beam along +x at the molecule-inlet port face, unless a
# thermal-* condition overrides it
SPAWN_X=0.03048; VX=24843.71; TSPAWN=0
case "$COND" in
  edep) A="--sigma-model edep" ;;   # 1/v-with-floor model (default before 55b5f6b)
  table-central|table-low|table-high)   # computed BaF+-Ne table (tracer default since 55b5f6b)
        A="--sigma-model table --sigma-table-variant ${COND#table-}" ;;
  const2e-18) A="--sigma-energy-dependent 0 --sigma 2.0e-18" ;;
  const1.5e-19) A="--sigma-energy-dependent 0 --sigma 1.5e-19" ;;
  thermal-const1.93e-18|thermal-const1.93e-18-x42|thermal-const1.93e-18-mid)
        # test: BaF+ already at the Ne temperature (23 K Maxwellian, no beam),
        # constant 1.93e-18 m^2, on the port axis. Start point:
        #   (no suffix) at the inlet channel's inner end (x=31.8mm): swept back
        #       out by the neon leaving through the inlet port
        #   -x42  10mm into the chamber (x=41.8mm), ~where 500 eV ions stop
        #   -mid  chamber centre (x=54.5mm), like molecules born in the cell
        A="--sigma-model const --sigma 1.93e-18"; VX=0; TSPAWN=23
        case "$COND" in *-x42) SPAWN_X=0.0418 ;; *-mid) SPAWN_X=0.0545 ;; *) SPAWN_X=0.0318 ;; esac ;;
  *) echo "unknown COND $COND"; exit 2 ;;
esac
RUNDIR="$HOME/orcd/scratch/cbgb-sparta-3d/fluor-scan/sccm$SCCM"
[ "$SCCM" = 4 ] && RUNDIR="$HOME/orcd/scratch/cbgb-sparta-3d/fluor-cell-3d-ne-sealed"
cd "$RUNDIR/post" || exit 2
[ -f frozen.grid ] || { echo "no frozen.grid -- run sbatch_fluor3d_scan_post.sh first"; exit 2; }
# retry a fast failure: concurrent jobs (or a node with a different CPU
# type) can collide while precompiling into the shared Julia depot
for attempt in 1 2 3; do
  t0=$(date +%s)
  "$SLURM_SUBMIT_DIR/.local/julia-1.9.4/bin/julia" --project="$SLURM_SUBMIT_DIR/tracer" \
    --threads="$SLURM_CPUS_PER_TASK" "$SLURM_SUBMIT_DIR/tracer/ParticleTracing3D.jl" \
    cell_fluor3d_sealed.surf frozen.grid -n "$NPART" -x "$SPAWN_X" -y 0.01905 -z 0.01905 \
    --vx "$VX" --vy 0 --vz 0 -T "$TSPAWN" -m 20.1797 -M 156.325 $A --saveall 1 --seed "$SEED" \
    > "tracer_${COND}_s$SEED.out.tmp" 2> "tracer_${COND}_s$SEED.err"
  rc=$?
  [ "$rc" -eq 0 ] || [ $(( $(date +%s) - t0 )) -gt 300 ] && break
  echo "attempt $attempt failed fast (rc $rc), retrying"; sleep $((30 * attempt))
done
# only a finished run becomes a .out, so a timed-out chunk is never pooled
[ "$rc" -eq 0 ] && mv "tracer_${COND}_s$SEED.out.tmp" "tracer_${COND}_s$SEED.out"
echo "SCCM=$SCCM COND=$COND SEED=$SEED NPART=$NPART julia exit $rc"
exit "$rc"
