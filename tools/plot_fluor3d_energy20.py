#!/usr/bin/env python3
"""keV injection at 20 SCCM with the extended computed cross-section table.

Two panels (no shared axis):
  left:  BaF+ extraction vs injection energy (log x), 95% Wilson intervals; all extracted
         ions and those leaving thermalized (< 0.1 eV), plus the previous edep model in gray
  right: stopping depth (first drop below 1 eV) past the port face vs energy, mean with the
         p10-p90 spread, from results/figures/fluor3d_20sccm_stopdepth.csv (40 ions each)

Writes OUT.png and OUT_summary.csv.
Usage: python tools/plot_fluor3d_energy20.py [--out results/figures/fluor3d_20sccm_energy_table]
"""
import argparse
import csv
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tracer_classify3d import classify, load  # noqa: E402
from plot_fluor3d_scan import INK, INK2, GRID, SURF, M_BAF, EV, E_THERMAL, wilson, style  # noqa: E402

POST = os.path.expanduser("~/orcd/scratch/cbgb-sparta-3d/fluor-scan/sccm20/post")
ENERGIES = [500, 750, 1000, 2000, 5000, 10000]
BLUE, GRAY = "#2a78d6", "#9a9994"   # reference palette slot 1; neutral reference
CENTRE_MM, OUTLET_MM = 54.5 - 30.48, 77.22 - 30.48


def tag(model, e):
    base = "table-central" if model == "table" else "edep"
    return base if e == 500 else f"{base}-e{e}"


def stats(model, e):
    t = tag(model, e)
    fs = glob.glob(f"{POST}/tracer_{t}.out") + glob.glob(f"{POST}/tracer_{t}_s*.out")
    if not fs:
        return None
    r = load(fs)
    tr, lost, _ = classify(r)
    E = 0.5 * M_BAF * (r[tr, 7:10] ** 2).sum(1) / EV
    p, lo, hi = wilson(tr.sum(), len(r))
    pt, lot, hit = wilson((E < E_THERMAL).sum(), len(r))
    return dict(model=model, energy_eV=e, n=len(r), extracted=int(tr.sum()), frac=p, lo=lo, hi=hi,
                frac_thermal=pt, lo_thermal=lot, hi_thermal=hit, back_out=float(lost.mean()),
                t_median_ms=float(np.median(r[tr, 11]) * 1e3) if tr.any() else np.nan)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/figures/fluor3d_20sccm_energy_table")
    a = ap.parse_args()
    rows = [s for m in ("table", "edep") for e in ENERGIES if (s := stats(m, e))]
    with open(a.out + "_summary.csv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    depth = list(csv.DictReader(open("results/figures/fluor3d_20sccm_stopdepth.csv")))

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 4.2), dpi=150)
    fig.patch.set_facecolor(SURF)
    for x in (ax, bx):
        style(x)
        x.set_xscale("log")
        x.set_xticks(ENERGIES)
        x.set_xticklabels(["0.5", "0.75", "1", "2", "5", "10"])
        x.minorticks_off()
        x.set_xlim(420, 12500)
        x.set_xlabel("Injection energy  [keV]", color=INK, fontsize=10)

    for model, col, lab in (("table", BLUE, "computed σ table"), ("edep", GRAY, "previous edep model (1/v, floored)")):
        d = [r for r in rows if r["model"] == model]
        e = np.array([r["energy_eV"] for r in d]); p = np.array([r["frac"] for r in d]) * 100
        err = [p - np.array([r["lo"] for r in d]) * 100, np.array([r["hi"] for r in d]) * 100 - p]
        if model == "table":
            ax.errorbar(e, p, yerr=err, color=col, lw=2, marker="o", ms=6, capsize=3, elinewidth=1,
                        mec=SURF, mew=1.5, label=f"{lab}: all extracted", zorder=3)
            pt = np.array([r["frac_thermal"] for r in d]) * 100
            errt = [pt - np.array([r["lo_thermal"] for r in d]) * 100, np.array([r["hi_thermal"] for r in d]) * 100 - pt]
            ax.errorbar(e, pt, yerr=errt, color=col, lw=1.5, ls="--", marker="o", ms=5, mfc=SURF, mec=col,
                        mew=1.2, capsize=2, elinewidth=0.8, label=f"{lab}: leaving thermalized (< {E_THERMAL:g} eV)", zorder=3)
        else:
            ax.plot(e, p, color=col, lw=1.5, ls="--", marker="o", ms=4, label=lab, zorder=2)
    ax.set_ylim(bottom=0)
    ax.set_ylabel("BaF⁺ extracted through outlet  [%]", color=INK, fontsize=10)
    ax.set_title("Extraction at 20 SCCM", color=INK, fontsize=10, loc="left")
    ax.legend(frameon=False, fontsize=7.5, labelcolor=INK2, loc="lower right")

    for model, col, lab in (("table", BLUE, "computed σ table"), ("edep", GRAY, "previous edep model")):
        d = [r for r in depth if r["model"] == model]
        e = np.array([float(r["energy_eV"]) for r in d]); m = np.array([float(r["depth_mean_mm"]) for r in d])
        lo = np.array([float(r["depth_p10_mm"]) for r in d]); hi = np.array([float(r["depth_p90_mm"]) for r in d])
        bx.errorbar(e, m, yerr=[m - lo, hi - m], color=col, lw=2 if model == "table" else 1.5,
                    ls="-" if model == "table" else "--", marker="o", ms=6 if model == "table" else 4,
                    capsize=3, elinewidth=1, mec=SURF, mew=1.5, label=lab, zorder=3)
    for y, txt in ((CENTRE_MM, "chamber centre"), (OUTLET_MM, "outlet aperture")):
        bx.axhline(y, color=INK2, lw=0.8, ls=":")
        bx.text(440, y + 0.8, txt, fontsize=7.5, color=INK2)
    bx.set_ylim(0, 50)
    bx.set_ylabel("Stopping depth past port face  [mm]", color=INK, fontsize=10)
    bx.set_title("Where ions first drop below 1 eV (mean, p10-p90)", color=INK, fontsize=10, loc="left")
    bx.legend(frameon=False, fontsize=7.5, labelcolor=INK2, loc="upper left", bbox_to_anchor=(0, 0.93))

    fig.suptitle("BaF⁺ injected along the port axis at 20 SCCM Ne (sealed cell)", color=INK, fontsize=11, x=0.01, ha="left")
    fig.text(0.01, 0.01, "Error bars left: 95% Wilson interval. Depths: 40 ions per energy (~±1 mm statistical). "
             "Electronic stopping and BaF⁺ dissociation are not modeled.", fontsize=7.5, color=INK2)
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    fig.savefig(a.out + ".png", facecolor=SURF)
    print(f"wrote {a.out}.png, {a.out}_summary.csv")


if __name__ == "__main__":
    main()
