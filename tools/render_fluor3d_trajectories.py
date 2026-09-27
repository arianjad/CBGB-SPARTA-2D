#!/usr/bin/env python3
"""Animate recorded BaF+ trajectories in the sealed fluor-cell-3d housing.

Input is a ParticleTracing3D.jl run made with --saveall 1 --trajprint K. Its
rows with a negative idx are the flight legs of ions 1..K; the positive rows
are every ion's final state. Ions spend microseconds slowing down from
500 eV, then milliseconds drifting and diffusing, so frames are spaced
LOGARITHMICALLY in time.

Each ion is a dot with a short trail, colored by its fate:
- extracted: left through the outlet
- wall hit: stuck
- back out the inlet
Ions that have left the box drop out of the frame; stuck ions stay where
they hit.

Writes:
  OUT.gif       the animation
  OUT_paths.png a still of all recorded paths
  OUT.npz       per-frame positions + fates, small enough to re-render locally
                (times, pos[ion, frame, xyz], alive[ion, frame], fate, walls path)

Usage: python tools/render_fluor3d_trajectories.py TRAJ.out --surf SURF --out OUT
         [--frames 160] [--tmin 1e-7] [--trail 6] [--xmax 0.11]
"""
import argparse
import os
import sys

import numpy as np
import pyvista as pv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tracer_classify3d import classify, load  # noqa: E402
from plot_field_3d_slices import read_surf  # noqa: E402
from export_vtk import build_wall_mesh  # noqa: E402

FATES = {  # reference palette slots 1-3, same order everywhere
    "extracted": "#2a78d6",
    "wall hit": "#eb6834",
    "back out the inlet": "#1baf7a",
}


def read_legs(path, k):
    """Leg rows (-idx ... xn yn zn ... time) for ions 1..k, streamed (file can be ~1 GB)."""
    chunks = []
    with open(path) as fh:
        buf = []
        for line in fh:
            if line.startswith("-"):
                buf.append(line)
                if len(buf) >= 500000:
                    chunks.append(np.loadtxt(buf, usecols=(0, 4, 5, 6, 11)))
                    buf = []
        if buf:
            chunks.append(np.loadtxt(buf, usecols=(0, 4, 5, 6, 11)))
    legs = np.vstack(chunks)
    legs[:, 0] = -legs[:, 0]
    return legs[legs[:, 0] <= k]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("traj")
    ap.add_argument("--surf", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--k", type=int, default=40)
    ap.add_argument("--frames", type=int, default=160)
    ap.add_argument("--tmin", type=float, default=1e-7)
    ap.add_argument("--trail", type=int, default=6)
    ap.add_argument("--xmax", type=float, default=0.11, help="crop view at this x (m); ions beyond are hidden")
    a = ap.parse_args()

    fin = load([a.traj])
    fin = fin[fin[:, 0] <= a.k]
    tr, lost, wall = classify(fin)
    fate = np.where(tr, 0, np.where(lost, 2, 1))          # 0 extracted, 1 wall hit, 2 back out inlet
    ids = fin[:, 0].astype(int)
    tend = fin[:, 11]
    legs = read_legs(a.traj, a.k)
    print(f"{len(legs)} legs for {len(ids)} ions; fates: "
          f"{(fate==0).sum()} extracted, {(fate==1).sum()} wall, {(fate==2).sum()} back out inlet")

    tmax = tend.max()
    times = np.logspace(np.log10(a.tmin), np.log10(tmax), a.frames)
    pos = np.full((len(ids), a.frames, 3), np.nan)
    alive = np.zeros((len(ids), a.frames), bool)
    for n, (i, te) in enumerate(zip(ids, tend)):
        L = legs[legs[:, 0] == i]
        t = np.concatenate([[0.0], L[:, 4], [te]])
        xyz = np.vstack([[0.03048, 0.01905, 0.01905], L[:, 1:4], fin[n, 4:7]])
        o = np.argsort(t, kind="stable"); t, xyz = t[o], xyz[o]
        j = np.clip(np.searchsorted(t, times, side="right") - 1, 0, len(t) - 1)
        pos[n] = xyz[j]
        # an extracted or back-out ion disappears once it has left; a stuck ion stays
        alive[n] = (times <= te) | (fate[n] == 1)
    visible = alive & (pos[:, :, 0] <= a.xmax)
    np.savez_compressed(a.out + ".npz", times=times, pos=pos, alive=visible, fate=fate, ids=ids,
                        fate_names=np.array(list(FATES)), surf=os.path.abspath(a.surf))

    points, tris, _ = read_surf(a.surf)
    walls = build_wall_mesh(points, tris)
    names = list(FATES)

    def fate_key(pl):
        # plain colored text lines: add_legend rescales its text to the box and clips it
        for k, nm in enumerate(names):
            pl.add_text(f"{nm} ({(fate == k).sum()})", position=(20, 20 + 26 * (len(names) - 1 - k)),
                        font_size=12, color=FATES[nm])

    def base_plotter():
        pl = pv.Plotter(off_screen=True, window_size=(1280, 720))
        pl.set_background("#fcfcfb")
        pl.add_mesh(walls, color="#b9b8b2", opacity=0.13, smooth_shading=True)
        pl.add_mesh(walls.extract_feature_edges(40), color="#8a8983", line_width=1)
        pl.camera_position = [(0.02, -0.17, 0.12), (0.058, 0.02, 0.016), (0, 0, 1)]
        pl.camera.zoom(1.6)
        return pl

    # still: every recorded path, downsampled for the render
    pl = base_plotter()
    for n in range(len(ids)):
        m = visible[n]
        if m.sum() > 1:
            pl.add_mesh(pv.lines_from_points(pos[n, m]), color=FATES[names[fate[n]]], line_width=2)
    fate_key(pl)
    pl.add_text("BaF+ paths, 17 SCCM Ne, 500 eV injection (log-time sampled)", font_size=11, color="#0b0b0b")
    pl.screenshot(a.out + "_paths.png")
    pl.close()

    # animation
    pl = base_plotter()
    pl.open_gif(a.out + ".gif", fps=12)
    fate_key(pl)
    dyn = []
    for f in range(a.frames):
        for act in dyn:
            pl.remove_actor(act)
        dyn = []
        for k, nm in enumerate(names):
            sel = np.where((fate == k) & visible[:, f])[0]
            if len(sel):
                dyn.append(pl.add_points(pos[sel, f], color=FATES[nm], point_size=11,
                                         render_points_as_spheres=True))
            for n in sel:
                lo = max(0, f - a.trail)
                seg = pos[n, lo:f + 1][visible[n, lo:f + 1]]
                if len(seg) > 1:
                    dyn.append(pl.add_mesh(pv.lines_from_points(seg), color=FATES[nm], line_width=2, opacity=0.6))
        t = times[f]
        label = f"t = {t*1e6:.2f} us" if t < 1e-3 else f"t = {t*1e3:.2f} ms"
        dyn.append(pl.add_text(f"{label}   (log time)   17 SCCM Ne, BaF+ injected at 500 eV; color = eventual fate",
                               font_size=11, color="#0b0b0b", position="upper_left"))
        pl.write_frame()
    pl.close()
    print(f"wrote {a.out}.gif, {a.out}_paths.png, {a.out}.npz")


if __name__ == "__main__":
    main()
