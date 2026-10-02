#!/usr/bin/env python3
"""Streaming version of field3d_freeze_avg.py for very large SPARTA grid dumps.

field3d_freeze_avg.py parses all wanted frames at once (and keeps their raw
text). For the 0.20 mm grid-convergence run (11.5M cells, ~2.8 GB/frame, 9
frames for 20 ms) that took >400 GB and >30 min. This version:
  - finds the frame offsets in the file's tail with `tail -c | grep -b`,
  - parses one frame at a time with pandas' C reader and keeps a running sum,
so memory is ~2 frames of float64 (a few GB).
Output is the same raw layout the tracer reads (SPARTA's 9-line ITEM: header,
then one row per cell). Empty duplicate frames that a resumed run re-dumps at
its start step are skipped, and the first non-empty copy of each step is kept.

Also prints the chamber summary line of tools/fluor3d_asymmetry.py (volume-
weighted <n>, T, <u,v,w> over the chamber box), computed on the averaged field.

Usage: python tools/field3d_freeze_stream.py FIELD_GRID --out OUT --span-steps S --max-frames K
"""
import argparse
import io
import os
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_field_3d_slices import NRHO_CEILING  # noqa: E402

NCOL = 15  # id xc yc zc xlo ylo zlo xhi yhi zhi nrho u v w T


def frame_offsets(path, k):
    """Byte offsets of the last frames' 'ITEM: TIMESTEP' lines (enough for k frames + margin)."""
    size = os.path.getsize(path)
    with open(path, "rb") as fh:
        head = fh.read(8 << 30 if size > 8 << 30 else size)
    second = head.find(b"ITEM: TIMESTEP", 1)
    frame_bytes = second if second > 0 else size
    del head
    want = min(size, int(frame_bytes * (2 * k + 4)))   # resumed jobs re-dump an empty duplicate per window
    out = subprocess.run(f"tail -c {want} '{path}' | grep -b -a '^ITEM: TIMESTEP'", shell=True,
                         capture_output=True, text=True, check=True).stdout
    base = size - want
    return [base + int(line.split(":", 1)[0]) for line in out.splitlines()], size


def read_frame(fh, off, end):
    """(timestep, header_text, DataFrame) for the frame in bytes [off, end)."""
    fh.seek(off)
    raw = fh.read(end - off)
    lines_end = 0
    for _ in range(9):
        lines_end = raw.index(b"\n", lines_end) + 1
    header = raw[:lines_end].decode()
    h = header.split("\n")
    step, ncell = int(h[1]), int(h[3])
    df = pd.read_csv(io.BytesIO(raw[lines_end:]), sep=" ", header=None, usecols=range(NCOL),
                     nrows=ncell, engine="c", dtype=np.float64)
    del raw
    if len(df) != ncell or df[0].isna().any():
        return step, None, None   # truncated or damaged frame
    return step, header, df


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("field_grid")
    ap.add_argument("--out", required=True)
    ap.add_argument("--span-steps", type=int, required=True)
    ap.add_argument("--max-frames", type=int, required=True)
    a = ap.parse_args()

    offs, size = frame_offsets(a.field_grid, a.max_frames)
    ends = offs[1:] + [size]
    # walk backwards from the end so the newest frames are read first, stopping once the span is covered
    total = ids = geom = header_out = None
    used, last = [], None
    with open(a.field_grid, "rb") as fh:
        for off, end in reversed(list(zip(offs, ends))):
            step, header, df = read_frame(fh, off, end)
            if df is None:
                print(f"skip truncated frame at step {step}")
                continue
            if last is None:
                last = step
            if step <= last - a.span_steps or len(used) >= a.max_frames:
                break
            if step in used or not np.any(df.values[:, 10:] != 0):
                continue   # duplicate or empty re-dump
            v = df.values
            d = v[np.argsort(v[:, 0], kind="stable")]
            if total is None:
                ids, geom, total, header_out = d[:, 0].copy(), d[:, 1:10].copy(), np.zeros((len(d), 5)), header
            elif not np.array_equal(d[:, 0], ids):
                sys.exit(f"cell ids differ at step {step}")
            total += d[:, 10:]
            used.append(step)
            print(f"frame {step}: {len(d)} cells", flush=True)
    avg = total / len(used)
    out = pd.DataFrame(np.column_stack([ids, geom, avg]))
    out[0] = out[0].astype(np.int64)
    with open(a.out, "w") as fo:
        fo.write(header_out)
        out.to_csv(fo, sep=" ", header=False, index=False, float_format="%.9g")
    print(f"averaged {len(used)} frames (steps {min(used)}-{max(used)}), {len(ids)} cells -> {a.out}")

    # chamber summary, same definition as tools/fluor3d_asymmetry.py
    xc, yc, zc = geom[:, 0]*1e3, geom[:, 1]*1e3, geom[:, 2]*1e3
    vol = np.prod(geom[:, 6:9] - geom[:, 3:6], axis=1) * 1e9
    n, u, v, w, T = (avg[:, i] for i in range(5))
    n = np.where(n < NRHO_CEILING, n, 0)
    N = n * vol
    inh = (xc > 30.48) & (xc < 82.55) & (yc > 0) & (yc < 44.45) & (zc > -6.35) & (zc < 38.1) & (n > 0)
    avgw = lambda q: (N*q)[inh].sum() / N[inh].sum()
    print(f"=== chamber: <n> {N[inh].sum()/vol[inh].sum():.3e} m^-3, T {avgw(T):.1f} K, "
          f"<u> {avgw(u):+.2f} <v> {avgw(v):+.2f} <w> {avgw(w):+.2f}")


if __name__ == "__main__":
    main()
