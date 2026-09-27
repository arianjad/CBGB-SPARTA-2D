#!/usr/bin/env python3
"""Pick 1-SCCM fill-in flows around the BaF+ extraction peak of the
fluor-cell-3d scan.

The peak is taken on the THERMALIZED extraction fraction (energy-dependent
cross section, < 0.1 eV at box exit), over flows >= 4 SCCM. At 1-2 SCCM the
"extracted" ions are hot ones that fly straight through, which isn't a beam.
Prints the integer flows strictly between the peak's two neighbouring
measured flows that haven't been run yet, capped at MAX nearest the peak.

Usage: python tools/fluor3d_pick_fine_flows.py [--max 8]
"""
import argparse
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tracer_classify3d import classify, load  # noqa: E402
from plot_fluor3d_scan import BASE, EXCLUDE, E_THERMAL, EV, M_BAF, post_dir  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--max", type=int, default=8)
a = ap.parse_args()
existing = sorted({4} | {int(d[4:]) for d in os.listdir(f"{BASE}/fluor-scan") if d.startswith("sccm")} - EXCLUDE)
res = {}
for f in existing:
    if f < 4:
        continue
    d = post_dir(f)
    fs = sorted(glob.glob(f"{d}/tracer_edep.out") + glob.glob(f"{d}/tracer_edep_s*.out"))
    if not fs:
        continue
    r = load(fs)
    tr, _, _ = classify(r)
    E = 0.5 * M_BAF * (r[tr, 7:10] ** 2).sum(1) / EV
    k = int((E < E_THERMAL).sum())
    res[f] = (k / len(r), k, len(r))
    print(f"# {f:4d} SCCM thermalized extraction {100*k/len(r):.2f}% ({k}/{len(r)})", file=sys.stderr)
flows = sorted(res)
peak = max(flows, key=lambda f: res[f][0])
i = flows.index(peak)
lo = flows[i - 1] if i > 0 else peak
hi = flows[i + 1] if i + 1 < len(flows) else peak + 5
cand = [f for f in range(lo + 1, hi) if f not in existing and f not in EXCLUDE]
cand = sorted(sorted(cand, key=lambda f: (abs(f - peak), f))[:a.max])
print(f"# peak {peak} SCCM, neighbours {lo} and {hi}", file=sys.stderr)
print(" ".join(map(str, cand)))
