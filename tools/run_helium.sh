#!/usr/bin/env bash
# Student B5 example: 4 SCCM, analytic prefill, 10 ms settling + 2 ms sampling.
# Optional environment: STEPS, RANKS, MDOT, FILLN, FNUM, DT, SEED.
set -euo pipefail
if [ "$#" -ne 1 ]; then echo "usage: $0 RUN_NAME" >&2; exit 2; fi
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
exec bash "$repo/cases/b5-lean/run_lean.sh" "$1" "${RANKS:-4}" in.he_b5_mflow \
  "MDOT=${MDOT:-1.1905173e-8},FILLN=${FILLN:-2.3438e21},FNUM=${FNUM:-2.5e17},DT=${DT:-1e-7},STEPS=${STEPS:-120000},WEIGHT=radius,SEED=${SEED:-8675309}" \
  plain b5-lean in.he_b5_mflow cell_b5.surf \
  tube_fwd.surf tube_rev.surf ap_fwd.surf ap_rev.surf
