#!/bin/bash
# Login-node driver for the fluor-cell-3d flow scan (5-30 SCCM extension).
# As each flow's neon run completes (restart.base.<TOTAL> exists), submit its
# field prep (sbatch_fluor3d_scan_post.sh, NPART=0), then its tracer jobs
# (sbatch_fluor3d_tracer_chunk.sh on mit_preemptable, 40 threads, 1 h), and
# replot once everything is finished (PREP_ONLY=1: field prep only). Idempotent via marker files in each
# post/ dir (.prep_submitted, .tracer_submitted: the job ids), so it can be
# killed and restarted at any point. Polls every 5 min; only runs squeue/ls.
#   nohup bash tools/wsl/fluor3d_scan_driver.sh > fluor3d_scan_driver.log 2>&1 &
cd "$HOME/CBGB-SPARTA-2D" || exit 1
B="$HOME/orcd/scratch/cbgb-sparta-3d/fluor-scan"
FLOWS="${FLOWS:-5 10 15 20 25 30}"
source tools/wsl/fluor3d_flow_params.sh
total() { flow_params "$1"; echo "$TOTAL"; }
# first-pass energy-dependent ions per flow (tracer cost per ion grows ~ density^2)
nedep() { local f=$1; if [ "$f" -le 5 ]; then echo 10000; elif [ "$f" -le 10 ]; then echo 5000;
          elif [ "$f" -le 15 ]; then echo 3000; elif [ "$f" -le 20 ]; then echo 2000;
          elif [ "$f" -le 25 ]; then echo 1500; else echo 1000; fi; }
TR=tools/wsl/sbatch_fluor3d_tracer_chunk.sh
# node17xx/node18xx crash inside LLVM while Julia precompiles CSV.jl there
# (a CPU-type problem, not a race, so retrying doesn't help)
# (a node range in --exclude must name only existing nodes, so list them)
BADNODES=$(sinfo -p mit_preemptable -h -N -o %N | grep -E "^node1[78][0-9][0-9]$" | sort -u | paste -sd,)
TOPT="-p mit_preemptable -t 01:00:00 --cpus-per-task=40 --mem=96G${BADNODES:+ --exclude=$BADNODES}"
while true; do
  alldone=1
  for F in $FLOWS; do
    P=$B/sccm$F/post; mkdir -p "$P"
    [ -f "$B/sccm$F/restart.base.$(total $F)" ] || { alldone=0; continue; }
    if [ ! -f "$P/.prep_submitted" ]; then
      j=$(sbatch --parsable --export=ALL,SCCM=$F,NPART=0 tools/wsl/sbatch_fluor3d_scan_post.sh) && echo "$j" > "$P/.prep_submitted"
      echo "$(date) $F SCCM prep $j"; alldone=0; continue
    fi
    if [ ! -f "$P/.tracer_submitted" ]; then
      prep=$(cat "$P/.prep_submitted")
      if squeue -j "$prep" -h 2>/dev/null | grep -q .; then alldone=0; continue; fi
      if [ ! -f "$P/frozen.grid" ]; then
        echo "$(date) $F SCCM prep $prep finished without frozen.grid -- resubmitting"
        rm -f "$P/.prep_submitted"; alldone=0; continue
      fi
      # PREP_ONLY=1: stop at the frozen field (tracers come from fluor3d_stats_driver.sh)
      if [ "${PREP_ONLY:-0}" = 1 ]; then
        : > "$P/.tracer_submitted"; echo "$(date) $F SCCM frozen field ready"; continue
      fi
      ids=$(sbatch --parsable $TOPT --export=ALL,SCCM=$F,SEED=201,NPART=$(nedep $F) $TR)
      for c in const2e-18 const1.5e-19; do
        ids="$ids $(sbatch --parsable $TOPT --export=ALL,SCCM=$F,SEED=42,NPART=1000,COND=$c $TR)"
      done
      echo "$ids" > "$P/.tracer_submitted"; echo "$(date) $F SCCM tracer $ids"; alldone=0; continue
    fi
    for j in $(cat "$P/.tracer_submitted"); do
      squeue -j "$j" -h 2>/dev/null | grep -q . && alldone=0
    done
  done
  [ "$alldone" = 1 ] && break
  sleep 300
done
[ "${PREP_ONLY:-0}" = 1 ] && exit 0
source .venv/bin/activate
python tools/plot_fluor3d_scan.py
column -s, -t results/figures/fluor3d_scan_summary.csv | cut -c1-140
