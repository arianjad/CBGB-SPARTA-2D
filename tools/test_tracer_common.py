#!/usr/bin/env python
"""Synthetic checks for the 2D first-crossing common reducer.

Run only under an allocated local test slot:
    conda run -n claude-code python tools/test_tracer_common.py

A terminal forward crossing plus hook-side forward and back crossings checks
the union and survival accounting on the projected-segment 2D path.
"""
import csv
import tempfile
from pathlib import Path

import tracer_common as tc


COLS = "idx x y z xnext ynext znext vx vy vz collides time".split()
RECS = "idx nlegs nfwd nback have cx cy vx vy vz tc nc mism".split()


def write_table(path, header, rows, sep):
    with path.open("w", newline="") as f:
        writer = csv.writer(f, delimiter=sep)
        writer.writerow(header)
        writer.writerows(rows)


def write_records(path, rows):
    write_table(path, RECS, rows, ",")


def write_legs(path, rows):
    write_table(path, COLS, rows, " ")


SURFS2D = """ITEM: TIMESTEP
0
ITEM: NUMBER OF SURFS
1
ITEM: BOX BOUNDS pp pp
0.0 0.10
0.0 0.02
ITEM: SURFS id v1x v1y v2x v2y
1 0.065 0.030 0.065 0.040
"""


def record(idx, have=0, nback=0):
    return [idx, 1, int(have), nback, int(have),
            0.001, 0.002, 1.0, 2.0, 10.0, 0.003, 1, 0.0]


def common_inputs(root, surfs, legs, records):
    root.mkdir()
    write_legs(root / "legs.out", legs)
    write_records(root / "recs.csv", records)
    return root, surfs


def test_2d_union(tmp):
    surfs = tmp / "cell2d.surfs"
    surfs.write_text(SURFS2D)
    # 1: terminal leg crosses forward.  2: hook already crossed and came back.
    # 3: hook crossed, then the terminal leg crosses back.  4: no crossing.
    legs = [
        [1, 0.000, 0.000, 0.040, 0.000, 0.000, 0.080, 0, 0, 10, 1, 0.001],
        [2, 0.000, 0.000, 0.070, 0.000, 0.000, 0.080, 0, 0, 10, 1, 0.004],
        [3, 0.000, 0.000, 0.080, 0.000, 0.000, 0.040, 0, 0, -10, 1, 0.004],
        [4, 0.000, 0.000, 0.080, 0.000, 0.000, 0.040, 0, 0, -10, 1, 0.004],
    ]
    root, surfs = common_inputs(tmp / "two", surfs, legs,
                                 [record(1), record(2, have=1, nback=1),
                                  record(3, have=1), record(4)])
    r = tc.score(root, surfs, 0.065, ref=root)
    assert (r["N_first"], r["N_back"], r["N_survive"]) == (3, 2, 1), r
    assert r["first_from_terminating_leg"] == 1, r
    assert r["first_from_hook"] == 2, r
    assert r["checks"]["endpoint_identity"], r
    assert r["reference_rows"]["identical"], r
    print("PASS 2D: terminal union and hook/terminal back-crossings")


def main():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        test_2d_union(tmp)
    print("1/1 PASS")


if __name__ == "__main__":
    main()
