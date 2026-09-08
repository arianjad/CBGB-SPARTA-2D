"""Checks for signed Cartesian trajectory projections and terminal contacts."""
import tempfile
from pathlib import Path

import numpy as np

import plot_trajectories as pt
from test_tracer_analyze import SURFS, write_out


with tempfile.TemporaryDirectory() as temp:
    root = Path(temp)
    surfs = root / "cell.surfs"
    surfs.write_text(SURFS)
    paths = root / "paths.out"
    # 1 crosses x=0 and reaches the plate.  Its raw final Cartesian point is
    # beyond the plate, so the visible terminal marker must use the same
    # (z,r)-clip fraction as tracer_analyze, interpolated on its x,y,z leg.
    # 2 leaves downstream with y < 0.  3 has no negative hook row because it
    # is ballistic (collides == 0); it is still a complete terminal-only path.
    # 4 has collisions but no hooks and must not acquire an invented track.
    write_out(paths, [
        (-1, -.004, .003, .020, .004, .003, .030, 0, 0, 100, 1, .0001),
        (1, .004, .003, .030, .006, .0045, .080, 0, 0, 100, 1, .0001),
        (-2, .001, -.002, .020, .001, -.002, .030, 0, 0, 100, 1, .0001),
        (2, .001, -.002, .030, .001, -.002, .120, 0, 0, 100, 1, .0001),
        (3, -.001, .002, .020, .001, -.002, .080, 0, 0, 100, 0, .0001),
        (4, .001, .001, .020, .002, .002, .030, 0, 0, 100, 2, .0001),
    ])

    tracks, fates, coll, *_ = pt.build(paths, surfs)

    # Preserve signed Cartesian positions: the z-x projection crosses the
    # centerline and the independent y coordinate does not collapse to r.
    assert tracks[1].shape[1] == 3
    assert tracks[1][0, 1] < 0 < tracks[1][1, 1], tracks[1]
    assert np.allclose(tracks[1][0], [.020, -.004, .003])
    assert np.allclose(tracks[2][0], [.020, .001, -.002])
    assert tracks[2][-1, 2] < 0, tracks[2]

    # Particle 1 reaches z=53 mm at f = (0.053 - 0.030)/(0.080 - 0.030).
    # The final marker must lie at that fraction of its Cartesian leg, rather
    # than on the radial-plane contact or at the raw free-flight endpoint.
    f = (.053 - .030) / (.080 - .030)
    want = np.array([.030, .004, .003]) + f * np.array([.050, .002, .0015])
    assert np.allclose(tracks[1][-1], want), (tracks[1][-1], want)
    assert fates[1] == "wall", fates
    assert fates[2] == "extracted", fates

    # A collision-free endpoint has a valid two-point path; a nonzero-collision
    # endpoint without negative hook rows has unknown history and is omitted.
    assert tracks[3].shape == (2, 3), tracks[3]
    assert np.allclose(tracks[3][-1], [.080, .001, -.002])
    assert fates[3] == "extracted" and coll[3] == 0
    assert 4 not in tracks and 4 not in fates

    # Existing output names remain valid.  No field.grid is present, so this
    # also covers the walls-and-paths-only branch.
    figures = root / "figures"
    pt.main([str(paths), str(surfs), "--outdir", str(figures)])
    for name in ("traj_overlay.png", "traj_aperture.png"):
        assert (figures / name).is_file() and (figures / name).stat().st_size > 5000
    extracted = root / "extracted"
    pt.main([str(paths), str(surfs), "--outdir", str(extracted), "--extracted-only"])
    for name in ("traj_extracted.png", "traj_extracted_aperture.png"):
        assert (extracted / name).is_file() and (extracted / name).stat().st_size > 5000

print("PASS signed Cartesian paths, clip-fraction terminals, and ballistic endpoints")
