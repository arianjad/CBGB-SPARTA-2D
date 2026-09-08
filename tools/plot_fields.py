"""Shared SPARTA wall and final-frame readers.

For helium figures, run tools/plot_fields_b5.py.
"""
import glob
import numpy as np

def wall_segments(rundir, pattern="cell*.surf"):
    """Polyline segments from a .surf in the run dir, or None.

    SPARTA .surf: 'Points' records are 'id x y'; 'Lines' records are
    'id p1 p2' or, when the deck reads with the type keyword,
    'id type p1 p2' -- the two point ids are always the last two fields.

    pattern selects the surface file; cell*.surf is the solid wall.
    Optional mesh surfaces are drawn separately from the wall loop.
    """
    paths = sorted(glob.glob("%s/%s" % (rundir, pattern)))
    if not paths:
        return None
    pts, segs, sect = {}, [], None
    for raw in open(paths[0]):
        tok = raw.split("#")[0].split()
        if not tok:
            continue
        if tok[0] in ("Points", "Lines"):
            sect = tok[0]
            continue
        if len(tok) == 2 and tok[1] in ("points", "lines"):
            continue
        if sect == "Points":
            pts[int(tok[0])] = (float(tok[1]), float(tok[2]))
        elif sect == "Lines":
            ids = [int(v) for v in tok]
            segs.append((pts[ids[-2]], pts[ids[-1]]))
    return segs or None

def last_frame(path):
    rows = None
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        if lines[i].startswith("ITEM: CELLS"):
            j = i + 1
            block = []
            while j < len(lines) and not lines[j].startswith("ITEM:"):
                block.append(lines[j].split())
                j += 1
            rows = block
            i = j
        else:
            i += 1
    return np.array(rows, dtype=float)
