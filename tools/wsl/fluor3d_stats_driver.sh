#!/bin/bash
# Login-node driver that tops up BaF+ tracer statistics for the fluor-cell-3d
# flow scan until each flow has >= TARGET extracted ions (condition COND,
# default edep; 150 is ~+-16% at 95%) or NMAX traced, then replots.
# Up to PERFLOW (default 2) chunks in flight per flow, MAXJOBS total.
# A chunk that times out halves that flow's chunk size (other failures don't).
# State in each post/ dir: .stats_jobs (id:seed per chunk in flight), .chunk,
# .seed, so it can be killed and restarted. Polls every 10 min; only
# squeue/sacct and a quick classify of the (small) tracer outputs.
#   nohup bash tools/wsl/fluor3d_stats_driver.sh > fluor3d_stats_driver.log 2>&1 &
cd "$HOME/CBGB-SPARTA-2D" || exit 1
source .venv/bin/activate
B="$HOME/orcd/scratch/cbgb-sparta-3d"
FLOWS="${FLOWS:-4 5 8 10 15 20 23 25 30}"
COND="${COND:-edep}"   # tracer condition to top up (see sbatch_fluor3d_tracer_chunk.sh)
# per-condition state files; edep keeps its original names
SFX=""; [ "$COND" != edep ] && SFX="_$COND"
TARGET="${TARGET:-150}"; NMAX="${NMAX:-40000}"; MAXJOBS="${MAXJOBS:-10}"
PERFLOW="${PERFLOW:-2}"   # chunks in flight per flow
MINCHUNK="${MINCHUNK:-100}"   # halving after timeouts stops here
# starting chunk size (ions per 1-h preemptable job; halved after a timeout),
# times CHUNKSCALE (e.g. 0.2 for 8-thread 15-min quicktest jobs via TOPT)
CHUNKSCALE="${CHUNKSCALE:-1}"
chunk0() { local f=$1 c; if [ "$f" -le 4 ]; then c=6000; elif [ "$f" -le 5 ]; then c=5000;
           elif [ "$f" -le 8 ]; then c=4000; elif [ "$f" -le 10 ]; then c=3000;
           elif [ "$f" -le 25 ]; then c=1200; else c=1000; fi
           python3 -c "print(max(100, int($c * $CHUNKSCALE)))"; }
TR=tools/wsl/sbatch_fluor3d_tracer_chunk.sh
# node17xx/node18xx crash in LLVM when Julia precompiles CSV.jl (CPU type, not
# a race, retrying doesn't help). --exclude ranges must name existing nodes only.
BADNODES=$(sinfo -p mit_preemptable -h -N -o %N | grep -E "^node1[78][0-9][0-9]$" | sort -u | paste -sd,)
# TOPT can be overridden, e.g. for quicktest when preemptable is congested:
#   TOPT="-p mit_quicktest -t 00:15:00 --cpus-per-task=8 --mem=32G" CHUNKSCALE=0.2
# Several option sets separated by ';' are used in turn, one per chunk (short
# walltimes on several partitions backfill sooner than one long request).
TOPT="${TOPT:--p mit_preemptable -t 01:00:00 --cpus-per-task=40 --mem=96G${BADNODES:+ --exclude=$BADNODES}}"
IFS=';' read -ra TOPTS <<<"$TOPT"; ti=0
RUNBASE="${RUNBASE:-fluor-scan}"   # e.g. fluor-scan-fine for the grid-convergence runs (exported to the chunks)
pdir() { [ "$1" = 4 ] && [ "$RUNBASE" = fluor-scan ] && echo "$B/fluor-cell-3d-ne-sealed/post" || echo "$B/$RUNBASE/sccm$1/post"; }
count() {  # prints "n extracted" for the pooled edep outputs of a post dir
  python - "$1" "$COND" <<'EOF'
import glob, sys
sys.path.insert(0, "tools")
from tracer_classify3d import load, classify
fs = sorted(glob.glob(f"{sys.argv[1]}/tracer_{sys.argv[2]}.out") + glob.glob(f"{sys.argv[1]}/tracer_{sys.argv[2]}_s*.out"))
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
    P=$(pdir $F); [ -f "$P/.stats_jobs$SFX" ] && inflight=$((inflight + $(wc -l < "$P/.stats_jobs$SFX")))
  done
  for F in $FLOWS; do
    P=$(pdir $F)
    [ -f "$P/frozen.grid" ] || { alldone=0; continue; }
    # wait for the first-pass tracer jobs (fluor3d_scan_driver.sh) to finish
    if [ -f "$P/.tracer_submitted" ]; then
      busy=0; for j in $(cat "$P/.tracer_submitted"); do squeue -j "$j" -h 2>/dev/null | grep -q . && busy=1; done
      [ "$busy" = 1 ] && { alldone=0; continue; }
    fi
    [ -f "$P/.chunk$SFX" ] || chunk0 $F > "$P/.chunk$SFX"
    [ -f "$P/.seed$SFX" ] || echo 300 > "$P/.seed$SFX"
    touch "$P/.stats_jobs$SFX"
    # reap finished chunks; halve the chunk size after a timeout
    keep=""
    while IFS=: read -r j s; do
      [ -z "$j" ] && continue
      if squeue -j "$j" -h 2>/dev/null | grep -q .; then keep="$keep$j:$s\n"; continue; fi
      if [ ! -f "$P/tracer_${COND}_s$s.out" ]; then
        st=$(sacct -j $j -n -X -o State | tr -d ' ')
        rm -f "$P/tracer_${COND}_s$s.out.tmp"
        if [ "$st" = TIMEOUT ]; then   # only a timeout means the chunk was too big
          c=$(( $(cat "$P/.chunk$SFX") / 2 )); [ "$c" -lt "$MINCHUNK" ] && c=$MINCHUNK; echo "$c" > "$P/.chunk$SFX"
        fi
        echo "$(date) $F SCCM chunk $j (seed $s) left no output ($st); chunk $(cat "$P/.chunk$SFX")"
      fi
    done < "$P/.stats_jobs$SFX"
    printf "%b" "$keep" > "$P/.stats_jobs$SFX"
    read n k <<<"$(count "$P")"
    if [ "$k" -ge "$TARGET" ] || [ "$n" -ge "$NMAX" ]; then
      [ -s "$P/.stats_jobs$SFX" ] && alldone=0; continue   # let in-flight chunks finish
    fi
    alldone=0
    while [ "$(wc -l < "$P/.stats_jobs$SFX")" -lt "$PERFLOW" ] && [ "$inflight" -lt "$MAXJOBS" ]; do
      s=$(cat "$P/.seed$SFX"); echo $((s + 1)) > "$P/.seed$SFX"
      c=$(cat "$P/.chunk$SFX")
      o=${TOPTS[$((ti % ${#TOPTS[@]}))]}; ti=$((ti + 1))
      j=$(sbatch --parsable $o --export=ALL,SCCM=$F,SEED=$s,NPART=$c,COND=$COND $TR) || break
      echo "$j:$s" >> "$P/.stats_jobs$SFX"; inflight=$((inflight + 1))
      echo "$(date) $F SCCM $COND: n=$n extracted=$k -> chunk $j seed $s NPART $c (${o%% -t*})"
    done
    [ -s "$P/.stats_jobs$SFX" ] && alldone=0
  done
  [ "$alldone" = 1 ] && break
  sleep 600
done
echo "$(date) all flows at target"
python tools/plot_fluor3d_scan.py
column -s, -t results/figures/fluor3d_scan_summary.csv | cut -c1-140
