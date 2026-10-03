#!/bin/bash
# Submit the plain-cylinder (flow x energy) scan as dependency-chained jobs, per flow:
#   neon DSMC (sbatch_cyl3d_neon.sh) -> freeze (sbatch_cyl3d_post.sh)
#   -> tracer jobs (sbatch_cyl3d_tracer.sh), each running EPERJOB energies.
# Nothing runs on the login node after submission. Run from the repo root:
#   FLOWS="5 10 15" ENERGIES="$(seq -s ' ' 200 200 5000)" bash tools/wsl/submit_cyl3d_campaign.sh
# Env: CONFIG (ports | ap5, default ports), FLOWS, ENERGIES, NPART (2000), SEED (1), EPERJOB (5),
#      RUNBASE (default per CONFIG: cyl-ports / cyl-scan),
#      NEON_OPT (extra sbatch options for the neon jobs; their tasks, memory and
#      walltime come from cyl3d_flow_params.sh per flow), TR_OPT (for the tracer jobs),
#      SKIP_NEON=1 (field already frozen: submit tracer jobs only).
set -euo pipefail
: "${FLOWS:?}"; : "${ENERGIES:?}"
NPART="${NPART:-2000}"; SEED="${SEED:-1}"; EPERJOB="${EPERJOB:-5}"; export CONFIG="${CONFIG:-ports}"
NEON_OPT="${NEON_OPT:-}"
source tools/wsl/cyl3d_flow_params.sh
cyl_flow_params 5; RUNBASE="${RUNBASE:-$RUNBASE_DEF}"
TR_OPT="${TR_OPT:---cpus-per-task=16 --mem=64G -t 1-00:00:00}"
read -ra EN <<<"$ENERGIES"
for F in $FLOWS; do
  dep=""
  if [ "${SKIP_NEON:-0}" != 1 ]; then
    cyl_flow_params "$F"
    j1=$(sbatch --parsable -n "$NP" --mem="${MEM_GB}G" -t "$TIME" $NEON_OPT -J "${CONFIG}-ne$F" --export=ALL,CONFIG=$CONFIG,SCCM=$F,RUNBASE=$RUNBASE tools/wsl/sbatch_cyl3d_neon.sh)
    j2=$(sbatch --parsable -d afterok:$j1 -J "${CONFIG}-post$F" --export=ALL,CONFIG=$CONFIG,SCCM=$F,RUNBASE=$RUNBASE tools/wsl/sbatch_cyl3d_post.sh)
    dep="-d afterok:$j2"
    echo "SCCM=$F neon $j1 ($NP tasks, ${MEM_GB}G, $TIME; cell $CELL_MM mm, DT $DT) post $j2"
  fi
  for ((i = 0; i < ${#EN[@]}; i += EPERJOB)); do
    E="${EN[*]:i:EPERJOB}"
    j=$(sbatch --parsable $dep $TR_OPT -J "${CONFIG}-tr$F" \
        --export=ALL,CONFIG=$CONFIG,SCCM=$F,RUNBASE=$RUNBASE,NPART=$NPART,SEED=$SEED,ENERGIES="$E" tools/wsl/sbatch_cyl3d_tracer.sh)
    echo "SCCM=$F tracer $j: $E eV"
  done
done
