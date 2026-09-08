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
tracer itself was given, so the figure cannot silently disagree with the run.
The greyscale flow background is optional and self-locating: the field
directory's provenance.json records the SPARTA run it was converted from, and
that run's field.grid supplies the per-cell |v| map.  With no run dir reachable
the walls and trajectories are drawn on their own.

Coordinates: tracer z is the SPARTA axial coordinate, and the tracer's
transverse pair (x, y) gives the SPARTA radius r = sqrt(x^2 + y^2).
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
from plot_fields import last_frame  # noqa: E402
from tracer_analyze import read_surfs, clip_legs, COLS  # noqa: E402

# Okabe-Ito, colourblind-safe; blue/vermillion match tracer_analyze.py.
C_EXTRACT, C_WALL, C_RADIAL = "#0072B2", "#D55E00", "#666666"
# Greys truncated at 62% so the fast plume never darkens past the trajectory
# hues; log norm because a cell body (~1 m/s) and its plume (~300 m/s) span
# two and a half decades.
CMAP = ListedColormap(plt.get_cmap("Greys")(np.linspace(0.0, 0.62, 256)))


def read_rows(outfile):
    """({pid: [(z, r), ...] leg starts}, {pid: endpoint row}).

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
                legs.setdefault(-idx, []).append((a[2], np.hypot(a[0], a[1])))
            else:
                ends[idx] = a
    if not legs:
        sys.exit("no --trajprint rows in %s: rerun the tracer with "
                 "--trajprint N (it is a particle count, not a flag)" % outfile)
    return legs, ends


def build(outfile, surfs):
    """Clipped polylines and fates, keyed by particle id."""
    seg, xap, bounds = read_surfs(surfs)
    legs, ends = read_rows(outfile)
    pids = sorted(p for p in legs if p in ends)
    missing = sorted(p for p in legs if p not in ends)
    if missing:
        print("WARNING: %d of %d traced particles have no endpoint row (%s) -- "
              "these died on a wall and the run lacked --saveall 1, so their "
              "final leg is unknown; they are dropped. Rerun with --saveall 1."
              % (len(missing), len(legs),
                 ",".join(str(p) for p in missing[:8])), file=sys.stderr)
    if not pids:
        sys.exit("no traced particle has an endpoint row: rerun with --saveall 1")

    # The terminating leg: start from the endpoint row, and its stored end is
    # the UNCLIPPED free-path endpoint (final radii of 250-300 mm show up on a
    # 19.9 mm domain), so it has to be clipped to the first crossing.
    p1 = np.array([[ends[p][2], np.hypot(ends[p][0], ends[p][1])] for p in pids])
    p2 = np.array([[ends[p][5], np.hypot(ends[p][3], ends[p][4])] for p in pids])
    hit, code = clip_legs(p1, p2, seg, bounds)

    tracks, fates, coll = {}, {}, {}
    for k, p in enumerate(pids):
        tracks[p] = np.vstack([np.array(legs[p]), p1[k], hit[k]])
        coll[p] = int(ends[p][9])
        if hit[k, 0] > xap:
            fates[p] = "extracted"
        elif code[k] == 1:
            fates[p] = "wall"
        else:
            fates[p] = "radial"
    return tracks, fates, coll, seg, xap, bounds


def background(ax, surfs):
    """Greyscale |v| map of the SPARTA run this field came from, if reachable."""
    prov = Path(surfs).parent / "provenance.json"
    if not prov.exists():
        return None
    try:
        rundir = json.loads(prov.read_text()).get("source_run_dir")
    except Exception:
        return None
    if not rundir or not os.path.exists(os.path.join(rundir, "field.grid")):
        return None
    data = last_frame(os.path.join(rundir, "field.grid"))
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


def draw(ax, tracks, colors, seg, surfs, xlim, ylim):
    pc = background(ax, surfs)
    # Walls come from the tracer's OWN geometry, mirrored about the axis, so
    # the drawing cannot disagree with what the particles collided against.
    for x1, y1, x2, y2 in seg:
        for s in (1, -1):
            ax.plot([x1, x2], [s * y1, s * y2], color="k", lw=1.4,
                    solid_capstyle="butt", zorder=3)
    for p, xy in tracks.items():
        c = colors[p]
        ax.plot(xy[:, 0], xy[:, 1], color=c, lw=1.0, alpha=0.55, zorder=4)
        ax.plot(xy[0, 0], xy[0, 1], "o", ms=5, mfc="white", mec=c, mew=1.2,
                zorder=6)
        ax.plot(xy[-1, 0], xy[-1, 1], "x", ms=6, color=c, mew=1.6, zorder=6)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.set_xlabel("axial position z (m)")
    ax.set_ylabel("radial position r (m)")
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
                    help="second panel window (default: auto, around the "
                         "extraction plane)")
    ap.add_argument("--extracted-only", action="store_true",
                    help="drop every trajectory that did not reach the "
                         "extraction plane and colour the survivors by their "
                         "BIRTH radius instead of by fate")
    a = ap.parse_args(argv)

    tracks, fates, coll, seg, xap, bounds = build(a.outfile, a.surfs)
    n_traced = len(tracks)
    if a.extracted_only:
        tracks = {p: xy for p, xy in tracks.items() if fates[p] == "extracted"}
        if not tracks:
            sys.exit("no traced particle was extracted: trace more particles "
                     "(--trajprint N with N >~ 20/extraction fraction)")
    handles, n_fate, n = legend_handles(fates)
    ncoll = [coll[p] for p in tracks]
    label = a.label or Path(a.outfile).stem
    if a.extracted_only:
        # tracks[p][0] is the first leg's START, i.e. the birthplace itself --
        # the same point the spawn marker is drawn at -- so no spawn.csv join.
        r0 = {p: xy[0, 1] for p, xy in tracks.items()}
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

    fig, ax = plt.subplots(figsize=(12, 3.6), constrained_layout=True)
    pc = draw(ax, tracks, colors, seg, a.surfs,
              (bounds[0, 0], bounds[0, 1]), (-rmax, rmax))
    if pc is not None:
        fig.colorbar(pc, ax=ax, label="He flow speed |v| (m/s)", pad=0.01,
                     fraction=0.03)
    if sm is None:
        ax.legend(handles=handles, loc="lower left", fontsize=8, ncol=2,
                  framealpha=0.9)
        title = "Tracer trajectories in the simulated He flow\n"
    else:
        fig.colorbar(sm, ax=ax, label="birth radius $r_0$ (m)", pad=0.01,
                     fraction=0.03)
        title = "Extracted tracer trajectories, coloured by birth radius\n"
    ax.set_title(title + sub, fontsize=10)
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
    fig, ax = plt.subplots(figsize=(8, 4.4), constrained_layout=True)
    pc = draw(ax, tracks, colors, seg, a.surfs, (zlo, zhi), (rlo, rhi))
    if pc is not None:
        fig.colorbar(pc, ax=ax, label="He flow speed |v| (m/s)", pad=0.01,
                     fraction=0.03)
    if sm is None:
        ax.legend(handles=handles[:n_fate], loc="lower right", fontsize=8,
                  framealpha=0.9)
    else:
        fig.colorbar(sm, ax=ax, label="birth radius $r_0$ (m)", pad=0.01,
                     fraction=0.03)
    ax.set_title("Extraction region (z = %.1f-%.1f mm)" % (zlo * 1e3, zhi * 1e3),
                 fontsize=10)
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
