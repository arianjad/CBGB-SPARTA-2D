#!/usr/bin/env python3
"""Overlay tracer trajectories on the helium flow field they were traced in.

Reads the per-leg rows the tracer writes under `--trajprint N` (negative-index
rows in its own stdout), so any ordinary run on any field, species and spawn
mode can be plotted -- there is no separate recorder to keep in sync.

Usage
-----
    python tools/plot_trajectories.py \\
        RUN_DIR/trace/paths.out RUN_DIR/field/cell.surfs \\
        [--outdir DIR] [--label TEXT] [--zoom ZLO ZHI RLO RHI]
        [--extracted-only]

`--extracted-only` keeps just the trajectories that reached the extraction
plane and colours each one by the radius it was BORN at, which is what makes
the birthplace-to-extraction correlation visible; the default fate colouring
(blue extracted / vermillion wall / grey radial) is unchanged.

Geometry, the domain box and the extraction plane all come from the .surfs the
tracer itself was given.  Geometry is shown only as a central meridional wall
outline; the figure does not perform a three-dimensional wall-intersection
calculation.  The greyscale flow background comes from the frozen field beside
the geometry.
For older outputs without that file, provenance selects its recorded source
timestep exactly. With neither source available, walls and paths are drawn alone.

Coordinates: tracer z is the SPARTA axial coordinate and the tracer records
both signed transverse coordinates (x, y).  The two panels show the actual
3D paths projected into the z-x and z-y planes.  Their helium backgrounds are
central meridional slices only, not gas sampled along an off-plane path.
"""
import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm, ListedColormap
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from field_io import FieldFormatError, read_complete_frames, select_frame  # noqa: E402
from tracer_analyze import read_surfs, clip_legs, COLS  # noqa: E402

# Okabe-Ito, colourblind-safe; blue/vermillion match tracer_analyze.py.
C_EXTRACT, C_WALL, C_RADIAL = "#0072B2", "#D55E00", "#666666"
# Greys truncated at 62% so the fast plume never darkens past the trajectory
# hues; log norm because a cell body (~1 m/s) and its plume (~300 m/s) span
# two and a half decades.
CMAP = ListedColormap(plt.get_cmap("Greys")(np.linspace(0.0, 0.62, 256)))


def read_rows(outfile):
    """({pid: [(z, x, y), ...] leg starts}, {pid: endpoint row}).

    `--trajprint` emits one row per leg in the ordinary 12-column format with a
    NEGATIVE index -- but the print sits in the `else` arm of propagate's loop,
    so the TERMINATING leg is returned and never printed.  The trajprint rows
    therefore stop where the final leg begins, and the death point exists only
    in the ordinary positive-index row.  Without `--saveall 1` the tracer emits
    those only for particles that reached the domain boundary, so wall deaths
    have no endpoint at all -- hence the check in build().
    """
    legs, ends = {}, {}
    with open(outfile) as f:
        for line in f:
            v = line.split()
            if len(v) != len(COLS):
                continue
            try:
                idx = int(v[0])
            except ValueError:
                continue
            a = [float(t) for t in v[1:]]
            if idx < 0:
                legs.setdefault(-idx, []).append((a[2], a[0], a[1]))
            else:
                ends[idx] = a
    return legs, ends


def build(outfile, surfs):
    """Clipped polylines and fates, keyed by particle id."""
    seg, xap, bounds = read_surfs(surfs)
    legs, ends = read_rows(outfile)
    # A particle that has no collisions reaches its terminal row without a
    # negative hook row.  It is a complete one-leg ballistic track.  A
    # nonzero-collision endpoint without hooks lacks its earlier path and is
    # deliberately not invented here.
    pids = sorted(set(p for p in legs if p in ends) |
                  {p for p, a in ends.items() if int(a[9]) == 0})
    missing = sorted(p for p in legs if p not in ends)
    if missing:
        print("WARNING: %d of %d traced particles have no endpoint row (%s) -- "
              "these died on a wall and the run lacked --saveall 1, so their "
              "final leg is unknown; they are dropped. Rerun with --saveall 1."
              % (len(missing), len(legs),
                 ",".join(str(p) for p in missing[:8])), file=sys.stderr)
    if not pids:
        sys.exit("no usable trajectory has an endpoint row: rerun the tracer "
                 "with --trajprint N and --saveall 1")

    # The terminating leg: start from the endpoint row, and its stored end is
    # the UNCLIPPED free-path endpoint (final radii of 250-300 mm show up on a
    # 19.9 mm domain), so it has to be clipped to the first crossing.
    p1xyz = np.array([[ends[p][2], ends[p][0], ends[p][1]] for p in pids])
    p2xyz = np.array([[ends[p][5], ends[p][3], ends[p][4]] for p in pids])
    p1 = np.column_stack([p1xyz[:, 0], np.hypot(p1xyz[:, 1], p1xyz[:, 2])])
    p2 = np.column_stack([p2xyz[:, 0], np.hypot(p2xyz[:, 1], p2xyz[:, 2])])
    hit, code = clip_legs(p1, p2, seg, bounds)

    # Fates intentionally remain the tracer analyzer's linear (z, r) model.
    # Use its contact fraction on the original Cartesian leg only for the
    # visible terminal marker; this is not a three-dimensional wall solve.
    d_zr = p2 - p1
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.sum((hit - p1) * d_zr, axis=1) / np.sum(d_zr * d_zr, axis=1)
    frac = np.clip(np.nan_to_num(frac, nan=1.0), 0.0, 1.0)
    terminal = p1xyz + frac[:, None] * (p2xyz - p1xyz)

    tracks, fates, coll = {}, {}, {}
    for k, p in enumerate(pids):
        tracks[p] = np.vstack([np.asarray(legs.get(p, ()), float).reshape(-1, 3),
                               p1xyz[k], terminal[k]])
        coll[p] = int(ends[p][9])
        if hit[k, 0] > xap:
            fates[p] = "extracted"
        elif code[k] == 1:
            fates[p] = "wall"
        else:
            fates[p] = "radial"
    return tracks, fates, coll, seg, xap, bounds


def background_data(surfs):
    """Return frozen data, or the exact source timestep named in provenance."""
    prov = Path(surfs).parent / "provenance.json"
    frozen = Path(surfs).parent / "field.grid"
    if frozen.is_file():
        return read_complete_frames(frozen).frames[-1].data
    if not prov.is_file():
        return None
    try:
        source = json.loads(prov.read_text())
        rundir = source.get("source_run_dir")
        step = source["source_timestep"]
    except (OSError, ValueError, KeyError):
        return None
    if not rundir or not os.path.exists(os.path.join(rundir, "field.grid")):
        return None
    try:
        data = select_frame(read_complete_frames(os.path.join(rundir, "field.grid")).frames,
                            timestep=step).data
    except FieldFormatError:
        return None
    return data


def background(ax, surfs):
    """Add the selected frozen |v| map when an extended grid is available."""
    data = background_data(surfs)
    if data is None:
        return None
    if data.shape[1] < 11:
        return None                      # legacy 7-column dump: no cell extents
    xlo, ylo, xhi, yhi = data[:, 3], data[:, 4], data[:, 5], data[:, 6]
    nrho, u, v = data[:, 7], data[:, 8], data[:, 9]
    good = nrho > 0
    speed = np.clip(np.hypot(u, v)[good], 1.0, None)
    verts = [[(a, s * b), (c, s * b), (c, s * d), (a, s * d)]
             for s in (1, -1)
             for a, b, c, d in zip(xlo[good], ylo[good], xhi[good], yhi[good])]
    pc = PolyCollection(verts, array=np.concatenate([speed, speed]),
                        cmap=CMAP, edgecolors="none",
                        norm=LogNorm(vmin=1.0, vmax=speed.max()))
    ax.add_collection(pc)
    return pc


def draw(ax, tracks, colors, seg, surfs, xlim, ylim, transverse_index,
         slice_label):
    pc = background(ax, surfs)
    ax.axhline(0.0, color="#444444", lw=0.7, ls="--", alpha=0.45, zorder=2)
    # The mirrored outline is a valid axisymmetric envelope in either central
    # meridional slice.  It is not a rendered off-plane wall intersection.
    for x1, y1, x2, y2 in seg:
        for s in (1, -1):
            ax.plot([x1, x2], [s * y1, s * y2], color="k", lw=1.4,
                    solid_capstyle="butt", zorder=3)
    for p, xyz in tracks.items():
        c = colors[p]
        ax.plot(xyz[:, 0], xyz[:, transverse_index], color=c, lw=1.0,
                alpha=0.55, zorder=4)
        ax.plot(xyz[0, 0], xyz[0, transverse_index], "o", ms=5, mfc="white", mec=c, mew=1.2,
                zorder=6)
        ax.plot(xyz[-1, 0], xyz[-1, transverse_index], "x", ms=6, color=c,
                mew=1.6, zorder=6)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.set_xlabel("axial position z (m)")
    ax.set_ylabel("signed transverse %s (m)" %
                  ("x" if transverse_index == 1 else "y"))
    ax.text(0.015, 0.985, slice_label, transform=ax.transAxes, va="top",
            ha="left", fontsize=7, color="#333333",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72,
                  "pad": 1.5})
    return pc


def legend_handles(fates):
    """One entry per fate that actually occurred, so the key never lies."""
    n = {k: sum(f == k for f in fates.values())
         for k in ("extracted", "wall", "radial")}
    label = {"extracted": "extracted past the geometry (%d)",
             "wall": "stuck on a wall (%d)",
             "radial": "left through the radial boundary (%d)"}
    h = [Line2D([], [], color=c, lw=1.4, label=label[k] % n[k])
         for k, c in (("extracted", C_EXTRACT), ("wall", C_WALL),
                      ("radial", C_RADIAL)) if n[k]]
    return h, len(h), n


def marker_handles():
    """Legend entries for markers whose meaning is independent of fate."""
    return [
        Line2D([], [], marker="o", ms=5, mfc="white", mec="#333333",
               mew=1.2, linestyle="none", label="birth"),
        Line2D([], [], marker="x", ms=6, color="#333333", mew=1.6,
               linestyle="none", label="estimated termination (z-r clipping)"),
    ]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outfile", help="tracer stdout from a --trajprint run")
    ap.add_argument("surfs", help="the .surfs the tracer was given")
    ap.add_argument("--outdir", required=True, help="directory for trajectory figures")
    ap.add_argument("--label", default=None,
                    help="subtitle text (default: the run's file name)")
    ap.add_argument("--zoom", nargs=4, type=float, default=None,
                    metavar=("ZLO", "ZHI", "RLO", "RHI"),
                    help="second-panel z and signed transverse-coordinate "
                         "window, applied to both projections (default: auto, "
                         "around the extraction plane)")
    ap.add_argument("--extracted-only", action="store_true",
                    help="drop every trajectory that did not reach the "
                         "extraction plane and colour the survivors by their "
                         "BIRTH radius instead of by fate")
    a = ap.parse_args(argv)

    tracks, fates, coll, seg, xap, bounds = build(a.outfile, a.surfs)
    n_traced = len(tracks)
    if a.extracted_only:
        tracks = {p: xyz for p, xyz in tracks.items() if fates[p] == "extracted"}
        if not tracks:
            sys.exit("no traced particle was extracted: trace more particles "
                     "(--trajprint N with N >~ 20/extraction fraction)")
    handles, _, n = legend_handles(fates)
    handles += marker_handles()
    ncoll = [coll[p] for p in tracks]
    label = a.label or Path(a.outfile).stem
    if a.extracted_only:
        # tracks[p][0] is the first leg's START, i.e. the birthplace itself --
        # the same point the spawn marker is drawn at -- so no spawn.csv join.
        r0 = {p: np.hypot(xyz[0, 1], xyz[0, 2]) for p, xyz in tracks.items()}
        rnorm = plt.Normalize(min(r0.values()), max(r0.values()))
        sm = plt.cm.ScalarMappable(norm=rnorm, cmap="viridis")
        colors = {p: sm.to_rgba(v) for p, v in r0.items()}
        sub = ("%s  |  %d extracted of %d traced, %d-%d collisions each  |  "
               "extraction plane z = %.1f mm"
               % (label, len(tracks), n_traced, min(ncoll), max(ncoll),
                  xap * 1e3))
    else:
        sm = None
        col = {"extracted": C_EXTRACT, "wall": C_WALL, "radial": C_RADIAL}
        colors = {p: col[fates[p]] for p in tracks}
        sub = ("%s  |  %d trajectories, %d-%d collisions each  |  extraction "
               "plane z = %.1f mm" % (label, len(tracks), min(ncoll),
                                      max(ncoll), xap * 1e3))

    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rmax = bounds[1, 1]
    written = []

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6), constrained_layout=True,
                             sharex=True, sharey=True)
    pc = draw(axes[0], tracks, colors, seg, a.surfs,
              (bounds[0, 0], bounds[0, 1]), (-rmax, rmax), 1,
              "He slice y = 0\nwall: meridional outline")
    draw(axes[1], tracks, colors, seg, a.surfs,
         (bounds[0, 0], bounds[0, 1]), (-rmax, rmax), 2,
         "He slice x = 0\nwall: meridional outline")
    if pc is not None:
        fig.colorbar(pc, ax=axes, label="He flow speed |v| (m/s)", pad=0.01,
                     fraction=0.03, shrink=0.60)
    if sm is None:
        axes[0].legend(handles=handles, loc="lower left", fontsize=8, ncol=1,
                       framealpha=0.9)
    else:
        fig.colorbar(sm, ax=axes, label="birth radius $r_0$ (m)", pad=0.01,
                     fraction=0.03, shrink=0.60)
        axes[0].legend(handles=handles[-2:], loc="lower left", fontsize=8,
                       framealpha=0.9)
    title = "3D trajectories: signed Cartesian projections\n"
    fig.suptitle(title + sub, fontsize=10)
    p = outdir / ("traj_extracted.png" if sm is not None
                  else "traj_overlay.png")
    fig.savefig(p, dpi=170)
    plt.close(fig)
    written.append(p)

    # Second panel: the extraction region.  Derived from the geometry rather
    # than hardcoded, so it follows a longer cell or a moved aperture.
    if a.zoom:
        zlo, zhi, rlo, rhi = a.zoom
    else:
        span = bounds[0, 1] - bounds[0, 0]
        zlo, zhi = xap - 0.18 * span, min(bounds[0, 1], xap + 0.10 * span)
        rlo, rhi = -0.35 * rmax, 0.35 * rmax
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), constrained_layout=True,
                             sharex=True, sharey=True)
    pc = draw(axes[0], tracks, colors, seg, a.surfs, (zlo, zhi), (rlo, rhi),
              1, "He slice y = 0\nwall: meridional outline")
    draw(axes[1], tracks, colors, seg, a.surfs, (zlo, zhi), (rlo, rhi), 2,
         "He slice x = 0\nwall: meridional outline")
    if pc is not None:
        fig.colorbar(pc, ax=axes, label="He flow speed |v| (m/s)", pad=0.01,
                     fraction=0.03, shrink=0.60)
    if sm is None:
        axes[0].legend(handles=handles, loc="lower right", fontsize=8,
                       framealpha=0.9)
    else:
        fig.colorbar(sm, ax=axes, label="birth radius $r_0$ (m)", pad=0.01,
                     fraction=0.03, shrink=0.60)
        axes[0].legend(handles=handles[-2:], loc="lower right", fontsize=8,
                       framealpha=0.9)
    fig.suptitle("Extraction region (z = %.1f-%.1f mm)" %
                 (zlo * 1e3, zhi * 1e3), fontsize=10)
    p = outdir / ("traj_extracted_aperture.png" if sm is not None
                  else "traj_aperture.png")
    fig.savefig(p, dpi=170)
    plt.close(fig)
    written.append(p)

    for p in written:
        print("wrote %s" % p)
    print("fates: %s" % ", ".join("%s %d" % kv for kv in n.items()))


if __name__ == "__main__":
    main()
