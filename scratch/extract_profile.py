"""Dump on-axis cells near the aperture (xc,yc,nrho,u,t) from a converged
field.grid to a small CSV for x-position scans, so the multi-GB file is read
once. Fast tail-only reader (run is already known to be converged).

Usage: python3 extract_profile.py <sccm> <run_dir> <out_csv>
"""
import csv
import os
import sys

sys.path.insert(0, os.path.expanduser("~/CBGB-SPARTA-2D/tools"))
from plot_fields_b5 import average_tail
import numpy as np

sccm = sys.argv[1]
run_dir = sys.argv[2]
out_csv = sys.argv[3]

path = os.path.join(run_dir, "field.grid")
d = average_tail(path, 0.25, fast_tail=True)
x = d["xc"]; y = d["yc"]
mask = (x > 0.045) & (x < 0.060) & (y < 0.001)
idx = np.where(mask)[0]

with open(out_csv, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["sccm", "xc", "yc", "nrho", "u", "t"])
    for i in idx:
        w.writerow([sccm, d["xc"][i], d["yc"][i], d["nrho"][i], d["u"][i], d["t"][i]])

print(f"sccm={sccm}: wrote {len(idx)} near-aperture cells -> {out_csv}")
