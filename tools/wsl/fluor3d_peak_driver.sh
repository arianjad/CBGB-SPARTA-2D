#!/bin/bash
# Login-node driver for the 1-SCCM fill-in around the BaF+ extraction peak.
#   1. Waits for fluor3d_scan_driver.sh and fluor3d_stats_driver.sh to finish.
#   2. Picks the fill-in flows (tools/fluor3d_pick_fine_flows.py: integer
#      flows between the thermalized-extraction peak's two neighbours, max 8).
#   3. Starts their neon runs as preemptable run-to-end jobs.
#   4. Runs fluor3d_scan_driver.sh (prep + first tracer pass), then
#      fluor3d_stats_driver.sh (top-up to >=150 extracted, replot) on them.
# Preemptable only (requeued and resumed on preemption). Quicktest chains for
# 8 flows would fill quicktest's 8-job queue cap and block the prep jobs.
#   nohup bash tools/wsl/fluor3d_peak_driver.sh > fluor3d_peak_driver.log 2>&1 &
cd "$HOME/CBGB-SPARTA-2D" || exit 1
source .venv/bin/activate
echo "$(date) waiting for the running scan/stats drivers"
while pgrep -u "$USER" -f "fluor3d_(scan|stats)_driver.sh" >/dev/null; do sleep 600; done
python tools/plot_fluor3d_scan.py >/dev/null
NEW=$(python tools/fluor3d_pick_fine_flows.py --max "${MAX:-8}")
echo "$(date) fill-in flows: $NEW"
[ -z "$NEW" ] && { echo "nothing to add"; exit 0; }
for F in $NEW; do
  j=$(sbatch --parsable -p mit_preemptable -t 2-00:00:00 \
        --export=ALL,SCCM=$F,MAXCHUNKS=0,CHAIN=0 tools/wsl/sbatch_fluor3d_scan.sh)
  echo "$(date) $F SCCM neon run $j"
done
FLOWS="$NEW" bash tools/wsl/fluor3d_scan_driver.sh
FLOWS="$NEW" bash tools/wsl/fluor3d_stats_driver.sh
echo "$(date) fill-in done"
