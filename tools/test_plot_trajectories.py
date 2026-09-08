"""A free-flight endpoint beyond the plate must not relabel a wall hit."""
import tempfile
from pathlib import Path

import plot_trajectories as pt
from test_tracer_analyze import SURFS, write_out

with tempfile.TemporaryDirectory() as temp:
    root = Path(temp)
    surfs = root / 'cell.surfs'
    surfs.write_text(SURFS)
    paths = root / 'paths.out'
    # Particle 1 hits the plate at r=4 mm before its raw endpoint reaches
    # z=80 mm. Particle 2 passes through the aperture and escapes downstream.
    write_out(paths, [
        (-1, .004, 0, .020, .004, 0, .030, 0, 0, 100, 1, .0001),
        (1, .004, 0, .030, .004, 0, .080, 0, 0, 100, 1, .0001),
        (-2, .001, 0, .020, .001, 0, .030, 0, 0, 100, 1, .0001),
        (2, .001, 0, .030, .001, 0, .120, 0, 0, 100, 1, .0001),
    ])
    tracks, fates, *_ = pt.build(paths, surfs)
    assert fates == {1: 'wall', 2: 'extracted'}, fates
    assert tracks[1][-1, 0] < .0535
    assert tracks[2][-1, 0] > .0535
print('PASS trajectory fates use clipped contacts')
