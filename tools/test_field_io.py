"""Synthetic complete-frame, frozen-field, and bounded-plot checks."""
import json
import tempfile
from pathlib import Path

import numpy as np

import field2tracer as f2t
import plot_fields_b5 as pfb
import plot_trajectories as pt
from field_io import FieldFormatError, read_complete_frames, resolved_dt, select_frame, write_frame


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

        crlf = root / "crlf.grid"
        original = frame(20_000, ROWS_FLOW).replace("\n", "\r\n").encode()
        crlf.write_bytes(original)
        copied = root / "crlf-copy.grid"
        write_frame(read_complete_frames(crlf).frames[0], copied)
        assert copied.read_bytes() == original
        assert read_complete_frames(copied).frames[0].timestep == 20_000

        averaged = pfb.average_tail(grid, frac=0, until_ms=5.0)
        assert averaged["steps"] == [20_000]
        assert averaged["dt"] == 2e-7
        assert pfb.frame_label(averaged) == "frames ending at steps 20000-20000 (4-4 ms)"
        assert np.allclose(pt.background_data(out / "cell.surfs")[:, 0], [1, 2])
        assert pfb.average_tail(out / "field.grid", frac=0)["dt"] == 2e-7

        # Sparse windows contribute less to plotted flow/temperature summaries,
        # while number density remains a time average.
        pooled = pfb.average_tail(grid, frac=1)
        assert np.isclose(pooled["nrho"][0], 5e20)
        assert np.isclose(pooled["u"][0], 82.0)
        assert np.isclose(pooled["v"][0], 8.2)
        assert np.isclose(pooled["t"][0], 8.5)
        sparse_t = root / "sparse-temperature.grid"
        no_t = [row.copy() for row in ROWS_LATE]
        no_t[0][10] = 0
        sparse_t.write_text(frame(20_000, ROWS_FLOW) + frame(40_000, no_t))
        recovered = pfb.average_tail(sparse_t, frac=1)
        assert np.isclose(recovered["nrho"][0], 5e20)
        assert np.isclose(recovered["u"][0], 82.0)
        assert np.isclose(recovered["t"][0], 4.0)

        # New SPARTA dumps append flow area. A zero-density open cell must be
        # retained in the axis profile and density denominator; a solid cell
        # with zero flow area must not enter either.
        vgrid = root / "field-with-vol.grid"
        vr = [ROWS_ZERO[0] + [1e-4],
              ROWS_FLOW[1] + [1e-4],
              [3, .045, .005, .04, 0, .05, .01, 0, 0, 0, 0, 0]]
        vgrid.write_text(frame(20_000, vr).replace(
            "ITEM: CELLS " + COLUMNS, "ITEM: CELLS " + COLUMNS + " vol"))
        with_vol = pfb.average_tail(vgrid, frac=0)
        assert np.array_equal(pfb.open_cells(with_vol), [True, True, False])
        assert np.array_equal(pfb.axis_profile(with_vol)[1], [0, 2e20])
        assert np.isnan(pfb.axis_profile(with_vol)[2][0])
        assert np.isnan(pfb.axis_profile(with_vol)[3][0])
        mean, count = pfb.region_mean(with_vol, .02, .05, .01)
        assert count == 2 and np.isclose(mean, 1e20)
        pfb.fields_figure(with_vol, {}, "zero-density gas", root / "vol.png")
        assert (root / "vol.png").is_file()
        other = {**with_vol, "nrho": with_vol["nrho"].copy()}
        other["nrho"][:2] = [1e20, 1e20]
        pfb.compare_figure(with_vol, other, "zero", "sampled", {},
                           root / "ratio.png", xrad=.035)
        assert (root / "ratio.png").is_file()
        bad_vol = root / "bad-vol.grid"
        bad_vol.write_text(vgrid.read_text().replace(" temp vol", " temp extra"))
        try:
            pfb.average_tail(bad_vol, frac=0)
        except FieldFormatError:
            pass
        else:
            raise AssertionError("non-vol 12th column was accepted")

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
