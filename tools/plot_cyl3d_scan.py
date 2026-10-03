#!/usr/bin/env python3
"""BaF+ extraction over (flow, beam energy) for the plain-cylinder cell (cases/cyl-3d).

Reads every RUNBASE/sccm<F>/post/tracer_<COND>-e<E>_s*.out (sbatch_cyl3d_tracer.sh),
pooling seed chunks, and writes:
  OUT_summary.csv    flow, energy, n, extracted and thermalized (< 0.1 eV at the
                     box exit) counts and % with 95% Wilson intervals
  OUT_surface.png    3D surfaces: extraction % (z) over flow x energy, all extracted
                     ions and the thermalized ones
  OUT_heatmap.png    the same two grids as annotated heatmaps (exact values)

Usage: python tools/plot_cyl3d_scan.py [--runbase cyl-scan] [--cond table-central]
                                       [--out results/figures/cyl3d_scan]
"""
import argparse
import glob
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tracer_classify3d import CYL_BOX as BOX, CYL_OUTLET_X as OUTLET_X, classify, load  # noqa: E402

BASE = os.path.expanduser("~/orcd/scratch/cbgb-sparta-3d")
AMU, EV = 1.66054e-27, 1.602177e-19
THERMAL_EV = 0.1
# reference sequential ramp (blue 100 -> 700), magnitude = extraction %
SEQ = LinearSegmentedColormap.from_list("seq_blue", ["#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#104281", "#0d366b"])
INK, MUTED = "#1f1f1e", "#6b6a66"


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan, np.nan
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return 100 * p, 100 * (c - h), 100 * (c + h)


def collect(runbase, cond):
    rows = []
    for d in glob.glob(f"{BASE}/{runbase}/sccm*/post"):
        f = int(re.search(r"sccm(\d+)", d).group(1))
        byE = {}
        for p in glob.glob(f"{d}/tracer_{cond}-e*_s*.out"):
            E = int(re.search(rf"tracer_{re.escape(cond)}-e(\d+)_s", p).group(1))
            byE.setdefault(E, []).append(p)
        for E, paths in byE.items():
            r = load(paths)
            if len(r) == 0:
                continue
            tr, lost, wall = classify(r, OUTLET_X, BOX)
            ke = 0.5 * 156.325 * AMU * np.sum(r[:, 7:10] ** 2, axis=1) / EV
            th = tr & (ke < THERMAL_EV)
            rows.append((f, E, len(r), tr.sum(), th.sum(), lost.sum(), wall.sum()))
    rows.sort()
    return rows


def grids(rows, col):
    flows = sorted({r[0] for r in rows}); ens = sorted({r[1] for r in rows})
    Z = np.full((len(ens), len(flows)), np.nan)
    for r in rows:
        Z[ens.index(r[1]), flows.index(r[0])] = 100 * r[col] / r[2]
    return np.array(flows), np.array(ens), Z


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runbase", default="cyl-scan")
    ap.add_argument("--cond", default="table-central")
    ap.add_argument("--out", default="results/figures/cyl3d_scan")
    a = ap.parse_args()
    rows = collect(a.runbase, a.cond)
    if not rows:
        sys.exit("no tracer outputs found")
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out + "_summary.csv", "w") as fo:
        fo.write("sccm,energy_eV,n,extracted,extracted_pct,extracted_lo,extracted_hi,"
                 "thermalized,thermalized_pct,thermalized_lo,thermalized_hi,back_out,wall\n")
        for f, E, n, k, kt, lost, wall in rows:
            e, t = wilson(k, n), wilson(kt, n)
            fo.write(f"{f},{E},{n},{k},{e[0]:.2f},{e[1]:.2f},{e[2]:.2f},{kt},{t[0]:.2f},{t[1]:.2f},{t[2]:.2f},{lost},{wall}\n")
    print(f"{len(rows)} (flow, energy) points -> {a.out}_summary.csv")

    panels = [(3, "all extracted ions"), (4, f"thermalized (< {THERMAL_EV} eV at exit)")]
    zmax = max(np.nanmax(grids(rows, c)[2]) for c, _ in panels) or 1.0

    flows, ens, _ = grids(rows, 3)
    if len(flows) > 1 and len(ens) > 1:
        fig = plt.figure(figsize=(13, 5.6))
        for i, (col, title) in enumerate(panels, 1):
            F, E, Z = grids(rows, col)
            X, Y = np.meshgrid(F, E / 1000)
            ax = fig.add_subplot(1, 2, i, projection="3d")
            ax.plot_surface(X, Y, np.ma.masked_invalid(Z), cmap=SEQ, vmin=0, vmax=zmax,
                            edgecolor="white", linewidth=0.3, antialiased=True)
            ax.set_xlabel("Ne flow (SCCM)", color=INK); ax.set_ylabel("BaF⁺ energy (keV)", color=INK)
            ax.set_zlabel("extraction (%)", color=INK); ax.set_zlim(0, zmax * 1.05)
            ax.set_title(title, color=INK, fontsize=11)
            ax.view_init(elev=28, azim=-130)
        fig.suptitle(f"BaF⁺ extraction, plain cylinder (30.5 mm bore, 45.5 mm long, 5 mm apertures), "
                     f"σ {a.cond}", color=INK, fontsize=12)
        fig.tight_layout(); fig.savefig(a.out + "_surface.png", dpi=150); plt.close(fig)
        print(f"-> {a.out}_surface.png")

    fig, axs = plt.subplots(1, 2, figsize=(max(8, 1.0 * len(flows) * 2 + 3), max(4, 0.32 * len(ens) + 2)),
                            squeeze=False)
    for ax, (col, title) in zip(axs[0], panels):
        F, E, Z = grids(rows, col)
        im = ax.imshow(Z, origin="lower", aspect="auto", cmap=SEQ, vmin=0, vmax=zmax)
        ax.set_xticks(range(len(F)), [str(f) for f in F]); ax.set_yticks(range(len(E)), [str(e) for e in E])
        ax.set_xlabel("Ne flow (SCCM)", color=INK); ax.set_ylabel("BaF⁺ energy (eV)", color=INK)
        ax.set_title(title, color=INK, fontsize=11)
        for (j, i), z in np.ndenumerate(Z):
            if np.isfinite(z):
                ax.text(i, j, f"{z:.1f}", ha="center", va="center", fontsize=7,
                        color="white" if z > 0.55 * zmax else INK)
        fig.colorbar(im, ax=ax, label="extraction (%)")
    fig.tight_layout(); fig.savefig(a.out + "_heatmap.png", dpi=150); plt.close(fig)
    print(f"-> {a.out}_heatmap.png")


if __name__ == "__main__":
    main()
