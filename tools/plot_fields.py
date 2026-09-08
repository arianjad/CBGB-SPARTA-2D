"""Shared SPARTA wall and complete-frame readers.

For helium figures, run tools/plot_fields_b5.py.
"""
import glob
import numpy as np

from field_io import read_complete_frames, select_frame

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

def last_frame(path, timestep=None):
    """Return the latest complete frame, or one exact complete timestep."""
    frames = read_complete_frames(path).frames
    return select_frame(frames, timestep=timestep).data
