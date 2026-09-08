#!/usr/bin/env python3
"""Convert a 2D SPARTA helium field and solid walls to molecule-tracer inputs.

Usage: python tools/field2tracer.py RUN_DIR --out OUTPUT_DIR [SELECTION]
Writes a frozen selected `field.grid`, copied `cell*.surf` files, DS2FF.DAT,
cell.surfs, and provenance.json. The default is the latest complete frame;
an explicit timestep must exist exactly.

DS2FF columns: axial x, radial y, T, number density, mass density,
axial velocity, radial velocity, and zero azimuthal velocity. The tracer's
Cartesian axial coordinate is z. Only cell*.surf supplies solid walls;
transparent ap_*/tube_* diagnostic surfaces must not enter tracer geometry.
Closing segments outside the domain are retained and clipped by domain tests.
Synthetic fields have >=100 cells and slight T/velocity variation to avoid
degenerate lookup-table ranges; they are software test inputs.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from field_io import (FieldFormatError, manifest_info, read_complete_frames,
                      resolved_dt, select_frame, write_frame)  # noqa: E402
from plot_fields import wall_segments  # noqa: E402

# Mass-density column convention; the tracer uses number density and drops
# this column. The SPARTA species mass is specified in cases/data/he.species.
HE_MASS = 6.6464731e-27


def fmt(v):
    return repr(float(v))


def _header(step, n, bounds, kind, cols):
    """The 9 header lines both files share (data starts at line 10)."""
    return "\n".join([
        "ITEM: TIMESTEP", str(step),
        "ITEM: NUMBER OF %s" % kind, str(n),
        "ITEM: BOX BOUNDS oo ao pp",
        "%s %s" % (fmt(bounds[0][0]), fmt(bounds[0][1])),
        "%s %s" % (fmt(bounds[1][0]), fmt(bounds[1][1])),
        "%s %s" % (fmt(bounds[2][0]), fmt(bounds[2][1])),
        "ITEM: %s" % cols,
    ]) + "\n"


def write_pair(outdir, step, bounds, cells, segs):
    """cells: iterable of (xc, yc, temp, nrho, u, v). segs: ((x1,y1),(x2,y2))."""
    os.makedirs(outdir, exist_ok=True)
    rows = ["%s %s %s %s %s %s %s 0.0" % (fmt(xc), fmt(yc), fmt(t), fmt(n),
                                          fmt(n * HE_MASS), fmt(u), fmt(v))
            for xc, yc, t, n, u, v in cells]
    with open(os.path.join(outdir, "DS2FF.DAT"), "w") as f:
        f.write(_header(step, len(rows), bounds, "CELLS",
                        "CELLS xc yc temp nrho massrho u v w"))
        f.write("\n".join(rows) + "\n")
    srows = ["%d %s %s %s %s" % (i + 1, fmt(a[0]), fmt(a[1]), fmt(b[0]), fmt(b[1]))
             for i, (a, b) in enumerate(segs)]
    with open(os.path.join(outdir, "cell.surfs"), "w") as f:
        f.write(_header(step, len(srows), bounds, "SURFS",
                        "SURFS id v1x v1y v2x v2y"))
        f.write("\n".join(srows) + "\n")
    return len(rows), len(srows)


def git_hash():
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                              text=True, cwd=os.path.dirname(os.path.abspath(__file__)),
                              check=True).stdout.strip()
    except Exception:
        return None


def _freeze_inputs(rundir, outdir, dumpname, frame):
    """Freeze the selected raw field and every cell surface beside conversion."""
    os.makedirs(outdir, exist_ok=True)
    write_frame(frame, os.path.join(outdir, dumpname))
    copied = []
    for source in sorted(Path(rundir).glob("cell*.surf")):
        target = Path(outdir) / source.name
        shutil.copy2(source, target)
        copied.append(target.name)
    return copied


def convert(rundir, outdir, dumpname="field.grid", timestep=None,
            until_step=None, until_ms=None, dt=None):
    path = os.path.join(rundir, dumpname)
    frames = read_complete_frames(path)
    source_status, manifest_dt = manifest_info(rundir)
    explicit_dt = resolved_dt(rundir, dt) if dt is not None else None
    selection_dt = (resolved_dt(rundir, dt) if until_ms is not None
                    else explicit_dt if explicit_dt is not None else manifest_dt)
    frame = select_frame(frames.frames, timestep=timestep, until_step=until_step,
                         until_ms=until_ms, dt=selection_dt)
    step, bounds, data = frame.timestep, frame.bounds, frame.data
    ncells, nframes = len(data), len(frames.frames)
    if data.shape[1] >= 11:
        xc, yc = data[:, 1], data[:, 2]
        nrho, u, v, temp = data[:, 7], data[:, 8], data[:, 9], data[:, 10]
    else:
        xc, yc = data[:, 1], data[:, 2]
        nrho, u, v, temp = data[:, 3], data[:, 4], data[:, 5], data[:, 6]
    segs = wall_segments(rundir)
    if not segs:
        sys.exit("no cell*.surf in %s: refusing to emit geometry-free input" % rundir)
    copied_surfaces = _freeze_inputs(rundir, outdir, dumpname, frame)
    nc, ns = write_pair(outdir, step, bounds, zip(xc, yc, temp, nrho, u, v), segs)

    run_id = None
    mpath = os.path.join(rundir, "manifest.json")
    if os.path.exists(mpath):
        try:
            run_id = json.load(open(mpath)).get("run_id")
        except Exception:
            pass
    prov = dict(source_run_dir=os.path.abspath(rundir), source_dump=dumpname,
                source_run_id=run_id, source_timestep=step,
                selection={"timestep": timestep, "until_step": until_step,
                           "until_ms": until_ms},
                source_observed_status=source_status,
                source_dt_s=selection_dt,
                header_number_of_cells=ncells, snapshots_in_dump=nframes,
                ignored_incomplete_tail=frames.ignored_incomplete_tail,
                frozen_dump=os.path.abspath(os.path.join(outdir, dumpname)),
                frozen_cell_surfs=copied_surfaces,
                emitted_cells=nc, emitted_wall_segments=ns,
                box_bounds=bounds, massrho_per_nrho=HE_MASS,
                converter="tools/field2tracer.py", converter_git_hash=git_hash(),
                date=datetime.now(timezone.utc).isoformat(), synthetic=False)
    json.dump(prov, open(os.path.join(outdir, "provenance.json"), "w"), indent=2)
    return prov


def synthetic(args, outdir):
    T, nrho, xlo, xhi, ylo, yhi, nx, ny = args
    nx, ny = int(nx), int(ny)
    if nx * ny < 100:
        sys.exit("--synthetic needs NX*NY >= 100 cells (tracer knn k=100); "
                 "use >= 200, e.g. 20 10")
    dx, dy = (xhi - xlo) / nx, (yhi - ylo) / ny
    cells = []
    for i in range(nx):
        x = xlo + (i + 0.5) * dx
        # near-uniform: a zero-width T range or zero flow makes the tracer's
        # lookup-table ranges degenerate (zero step -> ArgumentError/InexactError)
        t = T * (0.9975 + 0.005 * (i / (nx - 1) if nx > 1 else 0.5))
        for j in range(ny):
            cells.append((x, ylo + (j + 0.5) * dy, t, nrho, 0.1, 0.0))
    # one dummy wall strictly outside the domain: zero segments crashes the
    # tracer's maximum(geom[:,[2,4]]); this one can never be intersected.
    ox, oy = xhi + (xhi - xlo), yhi + 0.2 * (yhi - ylo)
    segs = [((xhi + 0.5 * (xhi - xlo), oy), (ox, oy))]
    bounds = [(xlo, xhi), (ylo, yhi), (-0.5, 0.5)]
    nc, ns = write_pair(outdir, 0, bounds, cells, segs)
    prov = dict(synthetic=True, T=T, nrho=nrho, box_bounds=bounds,
                nx=nx, ny=ny, flow_u=0.1, temperature_ramp=[T * 0.9975, T * 1.0025],
                emitted_cells=nc, emitted_wall_segments=ns,
                massrho_per_nrho=HE_MASS, converter="tools/field2tracer.py",
                converter_git_hash=git_hash(),
                date=datetime.now(timezone.utc).isoformat())
    json.dump(prov, open(os.path.join(outdir, "provenance.json"), "w"), indent=2)
    return prov


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("rundir", nargs="?")
    p.add_argument("--out", required=True)
    p.add_argument("--dump", default="field.grid")
    selection = p.add_mutually_exclusive_group()
    selection.add_argument("--timestep", type=int,
                           help="convert this exact complete timestep")
    selection.add_argument("--until-step", type=int,
                           help="use the latest complete frame at or before this step")
    selection.add_argument("--until-ms", type=float,
                           help="use the latest complete frame at or before this time (ms)")
    p.add_argument("--dt", type=float,
                   help="timestep size in seconds when --until-ms cannot read manifest.json")
    p.add_argument("--synthetic", nargs=8, type=float,
                   metavar=("T", "NRHO", "XLO", "XHI", "YLO", "YHI", "NX", "NY"))
    a = p.parse_args()
    try:
        prov = synthetic(a.synthetic, a.out) if a.synthetic else \
            convert(a.rundir or sys.exit("need a run dir or --synthetic"), a.out,
                    a.dump, a.timestep, a.until_step, a.until_ms, a.dt)
    except FieldFormatError as exc:
        sys.exit(str(exc))
    print("wrote %s: %d cells, %d wall segments" %
          (a.out, prov["emitted_cells"], prov["emitted_wall_segments"]))


if __name__ == "__main__":
    main()
