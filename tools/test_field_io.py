"""Synthetic complete-frame, frozen-field, and bounded-plot checks."""
import json
import tempfile
from pathlib import Path

import numpy as np

import field2tracer as f2t
import plot_fields_b5 as pfb
import plot_trajectories as pt
from field_io import FieldFormatError, read_complete_frames, resolved_dt, select_frame


COLUMNS = "id xc yc xlo ylo xhi yhi nrho u v temp"


def frame(step, rows, declared=None):
    count = len(rows) if declared is None else declared
    return ("ITEM: TIMESTEP\n%d\nITEM: NUMBER OF CELLS\n%d\n"
            "ITEM: BOX BOUNDS pp pp pp\n0 0.1\n0 0.02\n-0.5 0.5\n"
            "ITEM: CELLS %s\n%s" % (step, count, COLUMNS,
                                       "\n".join(" ".join(map(str, r)) for r in rows) + "\n"))


ROWS_ZERO = [[1, .025, .005, .02, 0, .03, .01, 0, 0, 0, 0],
             [2, .035, .005, .03, 0, .04, .01, 0, 0, 0, 0]]
ROWS_FLOW = [[1, .025, .005, .02, 0, .03, .01, 1e20, 10, 1, 4],
             [2, .035, .005, .03, 0, .04, .01, 2e20, 20, 2, 5]]
ROWS_LATE = [[1, .025, .005, .02, 0, .03, .01, 9e20, 90, 9, 9],
             [2, .035, .005, .03, 0, .04, .01, 8e20, 80, 8, 8]]


def write_surface(path):
    path.write_text("Points\n2 points\n1 0 0\n2 .1 0\nLines\n1 lines\n1 1 2\n")


def main():
    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        run = root / "run"
        run.mkdir()
        grid = run / "field.grid"
        # The last frame declares two rows but has only one: it is still being
        # appended and must not replace the complete 20,000-step field.
        grid.write_text(frame(0, ROWS_ZERO) + frame(20_000, ROWS_FLOW) +
                        frame(40_000, ROWS_LATE) + frame(60_000, ROWS_FLOW[:1], declared=2))
        write_surface(run / "cell_b5.surf")
        (run / "manifest.json").write_text(json.dumps({
            "status": "running", "runtime_controls": {"mpi_ranks": 1},
            "command": ["spa_mpi", "-in", "in.he", "-var", "DT", "2e-7"]}))

        frames = read_complete_frames(grid)
        assert [f.timestep for f in frames.frames] == [0, 20_000, 40_000]
        assert frames.ignored_incomplete_tail
        assert select_frame(frames.frames).timestep == 40_000
        assert select_frame(frames.frames, until_step=5_000).timestep == 0
        assert select_frame(frames.frames, until_ms=5.0, dt=2e-7).timestep == 20_000
        try:
            select_frame(frames.frames, timestep=60_000)
        except FieldFormatError:
            pass
        else:
            raise AssertionError("incomplete tail was selectable")

        out = root / "frozen"
        prov = f2t.convert(run, out, until_ms=5.0)
        frozen = read_complete_frames(out / "field.grid")
        assert [f.timestep for f in frozen.frames] == [20_000]
        assert (out / "cell_b5.surf").is_file()
        assert prov["source_observed_status"] == "running"
        assert prov["source_dt_s"] == 2e-7
        assert prov["ignored_incomplete_tail"]

        averaged = pfb.average_tail(grid, frac=0, until_ms=5.0)
        assert averaged["steps"] == [20_000]
        assert averaged["dt"] == 2e-7
        assert pfb.frame_label(averaged) == "frames ending at steps 20000-20000 (4-4 ms)"
        assert np.allclose(pt.background_data(out / "cell.surfs")[:, 0], [1, 2])
        assert pfb.average_tail(out / "field.grid", frac=0)["dt"] == 2e-7

        frozen_bytes = (out / "field.grid").read_bytes()
        # A live source advancing after conversion must not alter the frozen
        # background selected for trajectories.
        grid.write_text(frame(0, ROWS_ZERO) + frame(40_000, ROWS_LATE))
        assert (out / "field.grid").read_bytes() == frozen_bytes
        assert np.allclose(pt.background_data(out / "cell.surfs")[:, 7], [1e20, 2e20])

        # Every possible cut through the next valid frame is an in-progress
        # append, never a reason to reject the preceding complete frame.
        previous = frame(0, ROWS_ZERO)
        next_frame = frame(20_000, ROWS_FLOW)
        for cut in range(len(next_frame)):
            grid.write_text(previous + next_frame[:cut])
            parsed = read_complete_frames(grid)
            assert [f.timestep for f in parsed.frames] == [0], cut
            assert parsed.ignored_incomplete_tail or cut == 0

        grid.write_text(previous + "ITEM: not-a-frame")
        try:
            read_complete_frames(grid)
        except FieldFormatError:
            pass
        else:
            raise AssertionError("invalid ITEM line was mistaken for a partial header")

        unknown = root / "unknown"
        unknown.mkdir()
        try:
            resolved_dt(unknown)
        except FieldFormatError:
            pass
        else:
            raise AssertionError("unknown DT was accepted for time selection")

        malformed = root / "malformed.grid"
        malformed.write_text(frame(1, [[1, 0, 0, 0, 0, 0, 0, "bad", 0, 0, 0]]))
        try:
            read_complete_frames(malformed)
        except FieldFormatError:
            pass
        else:
            raise AssertionError("malformed complete frame was accepted")
    print("PASS complete-frame selection, freezing, and bounded plot inputs")


if __name__ == "__main__":
    main()
