#!/usr/bin/env python3
"""Usage: python tools/tracer_classify3d.py OUT [OUT ...]
Classify ParticleTracing3D.jl --saveall output (idx x y z xn yn zn vx vy vz collides time):
a box-face exit beyond the outlet counts as transmitted, a box exit at or
before it counts as lost, and an interior end counts as a wall hit."""
import sys

import numpy as np

OUTLET_X = 0.08255
BOX = (0.005, 0.200, -0.020, 0.070, -0.030, 0.080)
EPS = 1e-5


def load(paths):
    """Rows (idx x y z xn yn zn vx vy vz collides time) pooled from one or more outputs."""
    rows = []
    for path in paths:
        for line in open(path):
            p = line.split()
            if len(p) == 12 and not p[0].startswith("-"):
                try:
                    rows.append([float(q) for q in p])
                except ValueError:
                    pass
    return np.array(rows).reshape(-1, 12)


def classify(r):
    """Boolean masks (transmitted, lost, wall) for rows from load()."""
    xn, yn, zn = r[:, 4], r[:, 5], r[:, 6]
    at_box = np.zeros(len(r), bool)
    for c, lo, hi in ((xn, BOX[0], BOX[1]), (yn, BOX[2], BOX[3]), (zn, BOX[4], BOX[5])):
        at_box |= (abs(c - lo) < EPS) | (abs(c - hi) < EPS)
    return at_box & (xn > OUTLET_X), at_box & (xn <= OUTLET_X), ~at_box


def main():
    for path in sys.argv[1:]:
        r = load([path])
        tr, lost, stuck = classify(r)
        col, t = r[:, 10], r[:, 11]
        n = len(r); sp = np.linalg.norm(r[:, 7:10], axis=1)
        print(f"{path}: n={n}  transmitted {tr.sum()} ({100*tr.mean():.1f}%)  lost {lost.sum()} ({100*lost.mean():.1f}%)  "
              f"wall {stuck.sum()} ({100*stuck.mean():.1f}%)")
        print(f"  collides mean {col.mean():.0f} median {np.median(col):.0f} max {col.max():.0f}; ballistic {(col==0).sum()}; "
              f"exit speed mean {sp.mean():.0f} m/s; transit time median {np.median(t)*1e3:.3f} ms")
        if tr.any():
            print(f"  transmitted: exit speed mean {sp[tr].mean():.0f} m/s, collides mean {col[tr].mean():.0f}, "
                  f"time median {np.median(t[tr])*1e3:.3f} ms")
        if stuck.any():
            print(f"  wall-hit x positions (mm): quartiles {np.percentile(r[stuck, 4]*1e3, [5, 25, 50, 75, 95]).round(1)}")


if __name__ == "__main__":
    main()
