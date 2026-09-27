#!/bin/bash
# Login-node driver that tops up BaF+ tracer statistics for the fluor-cell-3d
# flow scan until each flow has >= TARGET extracted ions (energy-dependent
# condition; ~+-16% at 95%) or NMAX traced, then replots.
# - Keeps up to 2 chunks in flight per flow, and at most MAXJOBS in total.
# - A chunk that times out halves that flow's chunk size (other failures don't).
# - State lives in each post/ dir: .stats_jobs (id:seed per chunk in
#   flight), .chunk, .seed. So the script can be killed and restarted.
# - Polls every 10 min; runs only squeue/sacct and a small classify over
#   the (small) tracer outputs.
#   nohup bash tools/wsl/fluor3d_stats_driver.sh > fluor3d_stats_driver.log 2>&1 &
cd "$HOME/CBGB-SPARTA-2D" || exit 1
source .venv/bin/activate
B="$HOME/orcd/scratch/cbgb-sparta-3d"
FLOWS="${FLOWS:-4 5 8 10 15 20 23 25 30}"
TARGET="${TARGET:-150}"; NMAX="${NMAX:-40000}"; MAXJOBS="${MAXJOBS:-10}"
# starting chunk size (ions per 1-h preemptable job; halved after a timeout)
chunk0() { local f=$1; if [ "$f" -le 4 ]; then echo 6000; elif [ "$f" -le 5 ]; then echo 5000;
           elif [ "$f" -le 8 ]; then echo 4000; elif [ "$f" -le 10 ]; then echo 3000;
           elif [ "$f" -le 25 ]; then echo 1200; else echo 1000; fi; }
TR=tools/wsl/sbatch_fluor3d_tracer_chunk.sh
# node17xx/node18xx crash inside LLVM while Julia precompiles CSV.jl there
# (a CPU-type problem, not a race, so retrying doesn't help)
# (a node range in --exclude must name only existing nodes, so list them)
BADNODES=$(sinfo -p mit_preemptable -h -N -o %N | grep -E "^node1[78][0-9][0-9]$" | sort -u | paste -sd,)
TOPT="-p mit_preemptable -t 01:00:00 --cpus-per-task=40 --mem=96G${BADNODES:+ --exclude=$BADNODES}"
pdir() { [ "$1" = 4 ] && echo "$B/fluor-cell-3d-ne-sealed/post" || echo "$B/fluor-scan/sccm$1/post"; }
count() {  # prints "n extracted" for the pooled edep outputs of a post dir
  python - "$1" <<'EOF'
import glob, sys
sys.path.insert(0, "tools")
from tracer_classify3d import load, classify
fs = sorted(glob.glob(sys.argv[1] + "/tracer_edep.out") + glob.glob(sys.argv[1] + "/tracer_edep_s*.out"))
if not fs:
    print(0, 0)
else:
    r = load(fs); tr, _, _ = classify(r); print(len(r), int(tr.sum()))
EOF
}
while true; do
  alldone=1
  inflight=0
  for F in $FLOWS; do
    P=$(pdir $F); [ -f "$P/.stats_jobs" ] && inflight=$((inflight + $(wc -l < "$P/.stats_jobs")))
  done
  for F in $FLOWS; do
    P=$(pdir $F)
    [ -f "$P/frozen.grid" ] || { alldone=0; continue; }
    # wait for the first-pass tracer jobs (fluor3d_scan_driver.sh) to finish
    if [ -f "$P/.tracer_submitted" ]; then
      busy=0; for j in $(cat "$P/.tracer_submitted"); do squeue -j "$j" -h 2>/dev/null | grep -q . && busy=1; done
      [ "$busy" = 1 ] && { alldone=0; continue; }
    fi
    [ -f "$P/.chunk" ] || chunk0 $F > "$P/.chunk"
    [ -f "$P/.seed" ] || echo 300 > "$P/.seed"
    touch "$P/.stats_jobs"
    # reap finished chunks; halve the chunk size after a timeout
    keep=""
    while IFS=: read -r j s; do
      [ -z "$j" ] && continue
      if squeue -j "$j" -h 2>/dev/null | grep -q .; then keep="$keep$j:$s\n"; continue; fi
      if [ ! -f "$P/tracer_edep_s$s.out" ]; then
        st=$(sacct -j $j -n -X -o State | tr -d ' ')
        rm -f "$P/tracer_edep_s$s.out.tmp"
        if [ "$st" = TIMEOUT ]; then   # only a timeout means the chunk was too big
          c=$(( $(cat "$P/.chunk") / 2 )); [ "$c" -lt 100 ] && c=100; echo "$c" > "$P/.chunk"
        fi
        echo "$(date) $F SCCM chunk $j (seed $s) left no output ($st); chunk $(cat "$P/.chunk")"
      fi
    done < "$P/.stats_jobs"
    printf "%b" "$keep" > "$P/.stats_jobs"
    read n k <<<"$(count "$P")"
    if [ "$k" -ge "$TARGET" ] || [ "$n" -ge "$NMAX" ]; then
      [ -s "$P/.stats_jobs" ] && alldone=0; continue   # let in-flight chunks finish
    fi
    alldone=0
    while [ "$(wc -l < "$P/.stats_jobs")" -lt 2 ] && [ "$inflight" -lt "$MAXJOBS" ]; do
      s=$(cat "$P/.seed"); echo $((s + 1)) > "$P/.seed"
      c=$(cat "$P/.chunk")
      j=$(sbatch --parsable $TOPT --export=ALL,SCCM=$F,SEED=$s,NPART=$c $TR) || break
      echo "$j:$s" >> "$P/.stats_jobs"; inflight=$((inflight + 1))
      echo "$(date) $F SCCM: n=$n extracted=$k -> chunk $j seed $s NPART $c"
    done
    [ -s "$P/.stats_jobs" ] && alldone=0
  done
  [ "$alldone" = 1 ] && break
  sleep 600
done
echo "$(date) all flows at target"
python tools/plot_fluor3d_scan.py
column -s, -t results/figures/fluor3d_scan_summary.csv | cut -c1-140
