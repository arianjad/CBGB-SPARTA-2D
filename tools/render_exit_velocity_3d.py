#!/usr/bin/env python3
"""3D view of the exit velocities of the BaF+ ions that leave through the inlet (PyVista).

Same cell and camera as the trajectory GIF (tools/render_fluor3d_trajectories.py). Each
back-out ion is an arrow at the point where it leaves the simulation box, along its
velocity (length = velocity x ARROW_S), coloured by its transverse speed sqrt(v_y^2 + v_z^2).
Left: the GIF's perspective. Right: end-on, looking back into the inlet along +x, where
the arrows show only the transverse part (plume divergence shows as outward-pointing
arrows growing with radius).

Usage: python tools/render_exit_velocity_3d.py TAG --surf SURF --out OUT.png [--title ...] [--geom cylports]
(TAG = tracer output path without .out). Light: a single render, fine on a login node.
"""
import argparse
import os
import sys

import numpy as np
import pyvista as pv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from export_vtk import build_wall_mesh  # noqa: E402
from plot_field_3d_slices import read_surf  # noqa: E402
from render_fluor3d_trajectories import GEOMS  # noqa: E402
from tracer_classify3d import classify, load  # noqa: E402

ARROW_S = 2.0e-5          # s: a 100 m/s velocity draws as a 2 mm arrow


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tag")
    ap.add_argument("--surf", required=True)
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

    points, tris, _ = read_surf(a.surf)
    walls = build_wall_mesh(points, tris)
    cloud = pv.PolyData(p)
    cloud["v"] = v * ARROW_S
    cloud["transverse speed (m/s)"] = vperp
    clim = [0, float(np.percentile(vperp, 99))]
    vt = v.copy(); vt[:, 0] = 0.0                      # transverse part only, for the end-on view
    cloud_t = pv.PolyData(p)
    cloud_t["v"] = vt * ARROW_S
    cloud_t["transverse speed (m/s)"] = vperp

    pl = pv.Plotter(off_screen=True, shape=(1, 2), window_size=(1800, 720), border=False)
    for col, (cl, title) in enumerate([(cloud, "oblique, from outside the inlet: exit velocity"),
                                       (cloud_t, "end-on, looking back into the inlet: transverse part")]):
        pl.subplot(0, col)
        pl.set_background("#fcfcfb")
        pl.add_mesh(walls, color="#b9b8b2", opacity=0.13, smooth_shading=True)
        pl.add_mesh(walls.extract_feature_edges(40), color="#8a8983", line_width=1)
        glyphs = cl.glyph(orient="v", scale="v", factor=1.0, geom=pv.Arrow(tip_length=0.3, tip_radius=0.12, shaft_radius=0.04))
        pl.add_mesh(glyphs, scalars="transverse speed (m/s)", cmap="Blues", clim=clim, show_scalar_bar=(col == 1),
                    scalar_bar_args=dict(title="transverse speed (m/s)", color="#1f1f1e", fmt="%.0f",
                                         title_font_size=14, label_font_size=12,
                                         position_x=0.68, position_y=0.04, width=0.3, height=0.09))
        pl.add_text(title, font_size=11, color="#1f1f1e")
        if col == 0:
            # oblique, from outside the inlet: the exit disk in front, the cell behind it
            c = p.mean(axis=0)
            pl.camera_position = [(c[0] - 0.075, -0.085, 0.06), (c[0] + 0.02, 0.0, 0.0), (0, 0, 1)]
        else:
            c = p.mean(axis=0)
            pl.camera_position = [(c[0] - 0.12, 0.0, 0.0), (c[0], 0.0, 0.0), (0, 0, 1)]
            pl.camera.zoom(1.3)
    pl.subplot(0, 0)
    pl.add_text(f"{a.title}\n{back.sum()} ions out through the inlet; arrows at their exit points, 100 m/s = 2 mm",
                position="lower_left", font_size=10, color="#1f1f1e")
    pl.screenshot(a.out)
    pl.close()
    print("->", a.out)


if __name__ == "__main__":
    main()
