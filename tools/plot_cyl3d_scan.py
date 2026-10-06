#!/usr/bin/env python3
"""BaF+ extraction over (flow, beam energy) for the plain-cylinder cell (cases/cyl-3d).

Reads every RUNBASE/sccm<F>/post/tracer_<COND>-e<E>_s*.out (sbatch_cyl3d_tracer.sh),
pooling seed chunks, and writes:
  OUT_summary.csv    flow, energy, n, extracted and thermalized (< 0.1 eV at the
                     box exit) counts and % with 95% Wilson intervals
  OUT_surface.png    3D surfaces: extraction % (z) over flow x energy, all extracted
                     ions and the thermalized ones
  OUT_heatmap.png    the same two grids as annotated heatmaps (exact values)
  OUT_dedx.png       dE/dx per ion vs depth (the Bragg curve, from the tracer's
                     --hist energy deposition) at 1-5 keV, at the median finished flow
  OUT_stopdepth.png  mean depth where ions first drop below 1 eV (--ionlog) vs beam
                     energy, one line per flow

Usage: python tools/plot_cyl3d_scan.py [--config ports|ap5] [--runbase ...] [--cond table-central] [--out ...]
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
from tracer_classify3d import CYL_BOX as BOX, CYL_OUTLET_X, CYL_PORTS_OUTLET_X, classify, load  # noqa: E402

BASE = os.path.expanduser("~/orcd/scratch/cbgb-sparta-3d")
AMU, EV = 1.66054e-27, 1.602177e-19
THERMAL_EV = 0.1
# reference sequential ramp (blue 100 -> 700), magnitude = extraction %
SEQ = LinearSegmentedColormap.from_list("seq_blue", ["#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#104281", "#0d366b"])
INK, MUTED = "#1f1f1e", "#6b6a66"
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]   # reference categorical slots 1-5
SEQ_STEPS = ["#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95",
             "#104281", "#0d366b"]   # blue 250..700: ordinal by flow


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan, np.nan
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return 100 * p, 100 * (c - h), 100 * (c + h)


OUTLET_X = CYL_PORTS_OUTLET_X     # set from --config in main()


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
    ap.add_argument("--config", choices=["ports", "ap5", "ports-he"], default="ports")
    ap.add_argument("--runbase", default=None, help="default cyl-ports (ports) / cyl-scan (ap5)")
    ap.add_argument("--cond", default="table-central")
    ap.add_argument("--out", default=None, help="default results/figures/cyl3d_ports / cyl3d_scan")
    a = ap.parse_args()
    global OUTLET_X
    OUTLET_X = CYL_OUTLET_X if a.config == "ap5" else CYL_PORTS_OUTLET_X
    a.runbase = a.runbase or {"ports": "cyl-ports", "ap5": "cyl-scan", "ports-he": "cyl-ports-he"}[a.config]
    a.out = a.out or {"ports": "results/figures/cyl3d_ports", "ap5": "results/figures/cyl3d_scan",
                      "ports-he": "results/figures/cyl3d_ports-he"}[a.config]
    a.gas = "He" if a.config == "ports-he" else "Ne"
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

    panels = [(3, "all extracted ions", "extracted (% of injected)"),
              (4, f"extracted below {THERMAL_EV} eV", "extracted < 0.1 eV (% of injected)")]
    zmax = max(np.nanmax(grids(rows, c)[2]) for c, _, _ in panels) or 1.0

    flows, ens, _ = grids(rows, 3)
    # the surface uses only flows with every energy traced: a partly traced flow would be
    # drawn as straight "fins" across its missing energies (the heatmap shows partial flows)
    nE = len(ens)
    full = [r for r in rows if sum(1 for q in rows if q[0] == r[0]) == nE]
    sflows = sorted({r[0] for r in full})
    if len(sflows) > 1 and nE > 1:
        fig = plt.figure(figsize=(13, 5.6))
        for i, (col, title, zl) in enumerate(panels, 1):
            F, E, Z = grids(full, col)
            X, Y = np.meshgrid(F, E / 1000)
            ax = fig.add_subplot(1, 2, i, projection="3d")
            ax.plot_surface(X, Y, np.ma.masked_invalid(Z), cmap=SEQ, vmin=0, vmax=zmax,
                            edgecolor="white", linewidth=0.3, antialiased=True)
            ax.set_xlabel(f"{a.gas} flow (SCCM)", color=INK); ax.set_ylabel("BaF⁺ energy (keV)", color=INK)
            ax.set_zlabel(zl, color=INK); ax.set_zlim(0, zmax * 1.05)
            ax.set_title(title, color=INK, fontsize=11)
            ax.view_init(elev=28, azim=-130)
        ports = "5 mm apertures" if a.config == "ap5" else "old-cell ports"
        ports += ", 4 K helium" if a.gas == "He" else ", 23 K neon"
        fig.suptitle(f"BaF⁺ extraction, plain cylinder (30.5 mm bore, 45.5 mm long, {ports}), "
                     f"σ {a.cond}", color=INK, fontsize=12)
        fig.tight_layout(); fig.savefig(a.out + "_surface.png", dpi=150); plt.close(fig)
        print(f"-> {a.out}_surface.png")

    fig, axs = plt.subplots(1, 2, figsize=(max(8, 1.0 * len(flows) * 2 + 3), max(4, 0.32 * len(ens) + 2)),
                            squeeze=False)
    for ax, (col, title, zl) in zip(axs[0], panels):
        F, E, Z = grids(rows, col)
        im = ax.imshow(Z, origin="lower", aspect="auto", cmap=SEQ, vmin=0, vmax=zmax)
        ax.set_xticks(range(len(F)), [str(f) for f in F]); ax.set_yticks(range(len(E)), [str(e) for e in E])
        ax.set_xlabel(f"{a.gas} flow (SCCM)", color=INK); ax.set_ylabel("BaF⁺ energy (eV)", color=INK)
        ax.set_title(title, color=INK, fontsize=11)
        for (j, i), z in np.ndenumerate(Z):
            if np.isfinite(z):
                ax.text(i, j, f"{z:.1f}", ha="center", va="center", fontsize=7,
                        color="white" if z > 0.55 * zmax else INK)
        fig.colorbar(im, ax=ax, label=zl)
    fig.tight_layout(); fig.savefig(a.out + "_heatmap.png", dpi=150); plt.close(fig)
    print(f"-> {a.out}_heatmap.png")
    plot_diagnostics(a, sorted({r[0] for r in rows}))


def _tagfiles(runbase, f, cond, E, suffix):
    return sorted(glob.glob(f"{BASE}/{runbase}/sccm{f}/post/tracer_{cond}-e{E}_s*{suffix}"))


def plot_diagnostics(a, flows):
    """dE/dx (Bragg) curves and stopping depths from the --hist / --ionlog outputs."""
    import pandas as pd
    # dE/dx at the median finished flow, 1-5 keV
    have = [f for f in flows if _tagfiles(a.runbase, f, a.cond, 1000, "_xr.csv")]
    if have:
        f = have[len(have) // 2]
        fig, ax = plt.subplots(figsize=(7.5, 4.2))
        avail = sorted({int(re.search(r"-e(\d+)_s", q).group(1)) for q in
                        glob.glob(f"{BASE}/{a.runbase}/sccm{f}/post/tracer_{a.cond}-e*_s*_xr.csv")} & (set(range(1000, 5001, 1000)) if a.gas == "Ne" else set(range(500, 5001, 500))))
        for k, E in enumerate(avail[:: max(1, len(avail) // 5)][:5]):
            paths = _tagfiles(a.runbase, f, a.cond, E, "_xr.csv")
            if not paths:
                continue
            nions = 0; tot = None
            for pth in paths:
                n = int(re.search(r"over all (\d+) ions", open(pth).readline()).group(1)); nions += n
                d = pd.read_csv(pth, comment="#").groupby("x_m").edep_eV.sum()
                tot = d if tot is None else tot.add(d, fill_value=0)
            # rebin the 0.25 mm histogram to 1 mm: eV per mm per ion
            x0 = np.asarray(tot.index) * 1e3
            b = np.floor(x0).astype(int)
            xs = np.unique(b)
            y = np.array([tot.values[b == v].sum() for v in xs]) / nions / 1.0
            x = xs + 0.5
            ax.plot(x, y, color=CAT[k], lw=2, label=f"{E/1000:g} keV")
        ax.axvspan(-1.27, 0, color="#b9b8b2", alpha=0.35, lw=0); ax.axvline(45.47, color=MUTED, lw=1, ls=":")
        ax.set_xlim(-2, 52); ax.set_xlabel("depth x from the inlet wall (mm)", color=INK)
        ax.set_ylabel("dE/dx per ion (eV/mm, 1 mm bins)", color=INK)
        ax.legend(fontsize=8, frameon=False, title="beam energy", title_fontsize=8)
        ax.set_title(f"BaF⁺ energy loss vs depth, {f} SCCM {a.gas} (shaded: inlet wall; dotted: outlet wall, inner face)",
                     color=INK, fontsize=10)
        ax.grid(alpha=0.25)
        fig.tight_layout(); fig.savefig(a.out + "_dedx.png", dpi=150); plt.close(fig)
        print(f"-> {a.out}_dedx.png")
    # mean depth at first KE < 1 eV vs energy, per flow
    fig, ax = plt.subplots(figsize=(7.5, 4.2)); any_ = False
    for j, f in enumerate(flows):
        Es, xs, fr = [], [], []
        for pth in sorted(glob.glob(f"{BASE}/{a.runbase}/sccm{f}/post/tracer_{a.cond}-e*_s*_ions.csv")):
            E = int(re.search(r"-e(\d+)_s", pth).group(1))
            col = pd.read_csv(pth, comment="#").x_1eV
            x1 = col.dropna()
            if len(x1):
                Es.append(E); xs.append(x1.mean() * 1e3); fr.append(len(x1) / len(col))
        if Es:
            o = np.argsort(Es); any_ = True
            Es, xs, fr = np.array(Es)[o] / 1000, np.array(xs)[o], np.array(fr)[o]
            c = SEQ_STEPS[min(j, 9)]
            ok = fr >= 0.1     # an average over <10% of the ions is not representative: hollow, no line
            ax.plot(np.where(ok, Es, np.nan), np.where(ok, xs, np.nan), color=c, lw=2, marker="o", ms=3,
                    label=f"{f} SCCM")
            if (~ok).any():
                ax.plot(Es[~ok], xs[~ok], ls="none", marker="o", ms=4, mfc="none", mec=c)
    if any_:
        ax.set_xlabel("BaF⁺ beam energy (keV)", color=INK); ax.set_ylabel("mean depth of first KE < 1 eV (mm)", color=INK)
        ax.set_title("Stopping depth, averaged over the ions that drop below 1 eV (hollow: fewer than 10% do)",
                     color=INK, fontsize=10); ax.grid(alpha=0.25)
        ax.legend(fontsize=8, ncol=2, frameon=False)
        fig.tight_layout(); fig.savefig(a.out + "_stopdepth.png", dpi=150)
        print(f"-> {a.out}_stopdepth.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
