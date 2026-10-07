"""Dump on-axis downstream cells (xc,yc,nrho,u,t, aperture to near the domain
edge) from a converged field.grid to a small CSV, so downstream distances can
be scanned without re-reading the big file. Uses the fast tail-only reader.

Usage: python3 extract_downstream_profile.py <sccm> <run_dir> <out_csv>
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
mask = (x > 0.055) & (x < 0.150) & (y < 0.0015)
idx = np.where(mask)[0]

with open(out_csv, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["sccm", "xc", "yc", "nrho", "u", "t"])
    for i in idx:
        w.writerow([sccm, d["xc"][i], d["yc"][i], d["nrho"][i], d["u"][i], d["t"][i]])

print(f"sccm={sccm}: wrote {len(idx)} downstream cells -> {out_csv}")
