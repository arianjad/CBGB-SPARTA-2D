# Per-flow numerics for the fluor-cell-3d neon flow scan, shared by
# sbatch_fluor3d_scan.sh, sbatch_fluor3d_scan_post.sh and the login-node
# drivers. Source this file, then call:  flow_params <SCCM>
# It sets DT, DUMPFREQ (the averaging window, in steps), CHUNKS_QT (windows
# per 15-min quicktest job) and TOTAL (steps for 0.3 s simulated).
# Flows already run keep their exact historical values. Any other flow uses
# the same rule those were picked by:
#   DT = ~0.48 x the Ne-Ne mean collision time at the aperture-law chamber
#        density (2.6us x 4/SCCM), capped at 5e-7;
#   40000-step windows up to 10 SCCM, 60000 above.
flow_params() {
  local f=$1
  case "$f" in
    1|2|5|8)  DT=5e-7;   DUMPFREQ=40000; CHUNKS_QT=2 ;;
    4)        DT=5e-7;   DUMPFREQ=50000; CHUNKS_QT=2 ;;   # sealed reference run
    10)       DT=5e-7;   DUMPFREQ=40000; CHUNKS_QT=1 ;;
    15)       DT=3e-7;   DUMPFREQ=60000; CHUNKS_QT=1 ;;
    20)       DT=2.5e-7; DUMPFREQ=60000; CHUNKS_QT=1 ;;
    23|25)    DT=2e-7;   DUMPFREQ=60000; CHUNKS_QT=1 ;;
    30)       DT=1.6e-7; DUMPFREQ=60000; CHUNKS_QT=1 ;;
    58)       DT=1e-7;   DUMPFREQ=60000; CHUNKS_QT=1 ;;
    127)      DT=4e-8;   DUMPFREQ=60000; CHUNKS_QT=1 ;;
    *) read DT DUMPFREQ CHUNKS_QT <<<"$(python3 -c "
import math
f=float($f); dt=min(5e-7, 0.48*2.6e-6*4/f)
e=math.floor(math.log10(dt)); dt=math.floor(dt/10**(e-1))*10**(e-1)   # round down to 2 sig figs
print(f'{dt:.2g}', 40000 if f <= 10 else 60000, 2 if f <= 8 else 1)")" ;;
  esac
  TOTAL=$(python3 -c "print(int(round(0.3/$DT/$DUMPFREQ))*$DUMPFREQ)")
}
