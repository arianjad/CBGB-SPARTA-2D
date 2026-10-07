#!/usr/bin/env python3
"""Freeze a TIME-AVERAGED frame from a 3D SPARTA grid dump for
ParticleTracing3D.jl.

A single frame is already averaged over one dump window (fix ave/grid),
which is enough at 1-8 SCCM (20-25 ms windows) but not at high flow, where the
fluor-cell-3d scan uses 2.4-12 ms windows to fit mit_quicktest. This averages
the last K frames, matched by cell id, and writes one frame in the same raw
layout (9-line ITEM: header, then rows), so the tracer reads it unchanged.

The duplicate empty frame a resumed run re-dumps at its starting step is
dropped (the first copy of each timestep is kept).

Usage: python tools/field3d_freeze_avg.py FIELD_GRID --out OUT --span-steps S [--max-frames N]
  averages the frames within the last S steps (at least one frame);
  --max-frames reads only the last N frames of the file.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from field_io import read_complete_frames, read_tail_frames  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("field_grid")
    p.add_argument("--out", required=True)
    p.add_argument("--span-steps", type=int, required=True)
    p.add_argument("--max-frames", type=int, default=0,
                   help="read only the file's last N frames (enough to cover --span-steps); for multi-GB dumps")
    a = p.parse_args()

    if a.max_frames > 0:
        frames = list(read_tail_frames(a.field_grid, a.max_frames))
    else:
        seen, frames = set(), []
        for f in read_complete_frames(a.field_grid).frames:
            if f.timestep not in seen:
                seen.add(f.timestep)
                frames.append(f)
    last = frames[-1].timestep
    use = [f for f in frames if f.timestep > last - a.span_steps] or frames[-1:]
    ref = None
    total = None
    for f in use:
        d = f.data[np.argsort(f.data[:, 0], kind="stable")]
        if ref is None:
            ref, total = d, np.zeros_like(d)
        elif not np.array_equal(d[:, 0], ref[:, 0]):
            sys.exit(f"cell ids differ at step {f.timestep}")
        total += d
    avg = total / len(use)
    avg[:, 0] = ref[:, 0]
    header = "\n".join(use[-1].raw.split("\n")[:9])
    with open(a.out, "w") as out:
        out.write(header + "\n")
        np.savetxt(out, avg, fmt=["%d"] + ["%.9g"] * (avg.shape[1] - 1), delimiter=" ")
    print(f"averaged {len(use)} frames (steps {use[0].timestep}-{last}), {len(avg)} cells -> {a.out}")


if __name__ == "__main__":
    main()
