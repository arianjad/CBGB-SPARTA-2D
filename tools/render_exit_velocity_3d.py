#!/usr/bin/env python3
"""3D view of the transverse flow at the inlet exit (PyVista), one oblique angle.

Gas streamlines from the frozen DSMC field (frozen.grid: cell-centred u, v, w),
coloured by the gas's transverse speed sqrt(v^2 + w^2): seeded through the chamber
and traced both ways, so they show the gas converging on the inlet and fanning out
into the plume outside it. The BaF+ ions that leave through the inlet are drawn as
small markers at their exit points (where they leave the simulation box), coloured by
their own transverse speed on the same scale.

Usage: python tools/render_exit_velocity_3d.py TAG --surf SURF --field frozen.grid --out OUT.png [--title ...]
(TAG = tracer output path without .out). Light: one render of a coarse resampled field.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import pyvista as pv
from matplotlib import colormaps
from matplotlib.colors import ListedColormap

CMAP = ListedColormap(colormaps["Blues"](np.linspace(0.3, 1.0, 256)))   # light blue -> dark blue

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from export_vtk import build_wall_mesh  # noqa: E402
from plot_field_3d_slices import read_surf  # noqa: E402
from render_fluor3d_trajectories import GEOMS  # noqa: E402
from tracer_classify3d import classify, load  # noqa: E402

H = 4e-4                                   # m: resampling spacing for the velocity field


def velocity_grid(path, box):
    """Cell-centred (u, v, w) from a frozen.grid dump, resampled onto a uniform grid."""
    d = pd.read_csv(path, sep=r"\s+", skiprows=9, header=None, usecols=[1, 2, 3, 10, 11, 12, 13],
                    names=["x", "y", "z", "n", "u", "v", "w"])
    d = d[d.n > 0]
    cloud = pv.PolyData(d[["x", "y", "z"]].to_numpy())
    cloud["vel"] = d[["u", "v", "w"]].to_numpy()
    (x0, x1), (y0, y1), (z0, z1) = box
    dims = (int((x1 - x0) / H) + 1, int((y1 - y0) / H) + 1, int((z1 - z0) / H) + 1)
    g = pv.ImageData(dimensions=dims, spacing=(H, H, H), origin=(x0, y0, z0))
    g = g.interpolate(cloud, radius=2.0e-3, sharpness=2, strategy="null_value", null_value=0.0)
    g["transverse speed (m/s)"] = np.hypot(g["vel"][:, 1], g["vel"][:, 2])
    return g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tag")
    ap.add_argument("--surf", required=True)
    ap.add_argument("--field", required=True, help="frozen.grid of the run")
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default="")
    ap.add_argument("--geom", default="cylports", choices=list(GEOMS))
    a = ap.parse_args()
    G = GEOMS[a.geom]
    r = load([a.tag + ".out"])
    tr, back, wall = classify(r, *G["classify"])
    p = r[back, 4:7]
    v = r[back, 7:10]
    vperp = np.hypot(v[:, 1], v[:, 2])

    g = velocity_grid(a.field, ((-0.0224, 0.0470), (-0.0160, 0.0160), (-0.0160, 0.0160)))
    # seeds: discs in and just inside the inlet bore, where the flow is fast and well resolved.
    # Traced forward out into the plume, and only a short way back into the chamber, where the
    # mean flow (~0.2 m/s) is buried in DSMC noise and streamlines would just wander.
    def disc(x, rmax, nr, nth):
        rr, th = np.meshgrid(np.linspace(0.0, rmax, nr), np.linspace(0, 2 * np.pi, nth, endpoint=False))
        pts = np.c_[np.full(rr.size, x), (rr * np.cos(th)).ravel(), (rr * np.sin(th)).ravel()]
        return np.unique(np.round(pts, 7), axis=0)
    src = pv.PolyData(np.vstack([disc(-0.0006, 0.0011, 4, 12), disc(0.0015, 0.004, 4, 16)]))
    fwd = g.streamlines_from_source(src, vectors="vel", integration_direction="forward",    # along the flow: out
                                    max_length=0.035, initial_step_length=0.5, max_steps=20000, terminal_speed=0.5)
    bwd = g.streamlines_from_source(src, vectors="vel", integration_direction="backward",   # upstream: into the chamber
                                    max_length=0.008, initial_step_length=0.5, max_steps=20000, terminal_speed=0.5)
    lines = fwd.merge(bwd) if bwd.n_points else fwd
    clim = [0, float(np.percentile(np.r_[vperp, g["transverse speed (m/s)"]], 99.5))]

    points, tris, _ = read_surf(a.surf)
    walls = build_wall_mesh(points, tris)
    pl = pv.Plotter(off_screen=True, window_size=(1400, 800))
    pl.set_background("#fcfcfb")
    pl.add_mesh(walls, color="#b9b8b2", opacity=0.13, smooth_shading=True)
    pl.add_mesh(walls.extract_feature_edges(40), color="#8a8983", line_width=1)
    sb = dict(title="transverse speed (m/s)", color="#1f1f1e", fmt="%.0f", title_font_size=16,
              label_font_size=13, position_x=0.62, position_y=0.05, width=0.33, height=0.08)
    if lines.n_points:
        pl.add_mesh(lines.tube(radius=1.2e-4), scalars="transverse speed (m/s)", cmap=CMAP, clim=clim,
                    scalar_bar_args=sb)
    ions = pv.PolyData(p)
    ions["transverse speed (m/s)"] = vperp
    pl.add_mesh(ions, scalars="transverse speed (m/s)", cmap=CMAP, clim=clim, show_scalar_bar=not lines.n_points,
                scalar_bar_args=sb, render_points_as_spheres=True, point_size=7)
    c = p.mean(axis=0) if len(p) else np.array([-0.02, 0, 0])
    pl.camera_position = [(c[0] - 0.075, -0.085, 0.06), (c[0] + 0.02, 0.0, 0.0), (0, 0, 1)]
    pl.add_text(f"{a.title}: neon streamlines (lines) and BaF+ exit points (dots), colored by transverse speed",
                font_size=11, color="#1f1f1e")
    pl.add_text(f"{back.sum()} ions leave through the inlet; dots are where they leave the simulation box",
                position="lower_left", font_size=10, color="#1f1f1e")
    pl.screenshot(a.out)
    pl.close()
    print("->", a.out, f"({lines.n_lines} streamlines)")


if __name__ == "__main__":
    main()
