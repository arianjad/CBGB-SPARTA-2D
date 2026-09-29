#!/usr/bin/env python3
"""Plots for the fluor-cell-3d BaF+ flow-rate scan (sealed cell):

  OUT_extraction.png  -- extraction fraction vs flow (SCCM, log x) for the
                         three cross-section conditions, 95% Wilson intervals
  OUT_pumpout.png     -- one histogram per flow of the time transmitted ions
                         take from injection to leaving through the outlet
                         (energy-dependent condition)
  OUT_energy.png      -- one histogram per flow of extracted ions' kinetic
                         energy where they leave the box (x=0.2m, ~117mm past
                         the outlet), energy-dependent condition
  OUT_summary.csv     -- the numbers behind these

Reads every RUNDIR/post/tracer_<cond>*.out, pooling all seed chunks
(sbatch_fluor3d_scan_post.sh + sbatch_fluor3d_tracer_chunk.sh). Flows with
no tracer output yet are skipped.

Usage: python tools/plot_fluor3d_scan.py [--out results/figures/fluor3d_scan]
"""
import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tracer_classify3d import classify, load  # noqa: E402

BASE = os.path.expanduser("~/orcd/scratch/cbgb-sparta-3d")
EXCLUDE = {127}   # dropped for the 500 eV sets: ions don't leave through the outlet there
                  # (kept for the thermal set, as an upper bracket)
FLOWS_ALL = sorted({4} | {int(d[4:]) for d in os.listdir(f"{BASE}/fluor-scan") if d[4:].isdigit()})
FLOWS = [f for f in FLOWS_ALL if f not in EXCLUDE]
UNDER_RESOLVED = {f for f in FLOWS if f >= 10}   # Ne mean free path < half the 0.41mm chamber cell
CONDSETS = {  # (file tag, label, categorical slot colour -- reference palette slots 1-3); first = primary
    "edep": [
        ("edep", "energy-dependent σ", "#2a78d6"),
        ("const2e-18", "constant σ = 2.0e-18 m²", "#eb6834"),
        ("const1.5e-19", "constant σ = 1.5e-19 m²", "#1baf7a"),
    ],
    "table": [  # computed BaF+-Ne cross-section table (tracer default since 55b5f6b)
        ("table-central", "computed σ table, central", "#2a78d6"),
        ("table-low", "computed σ table, low variant", "#eb6834"),
        ("table-high", "computed σ table, high variant", "#1baf7a"),
    ],
}
CONDSETS["thermal"] = [  # test: BaF+ injected already at 23 K (no beam), constant 1.93e-18 m^2
    ("thermal-const1.93e-18-mid", "thermal start, chamber centre (x = 54.5 mm)", "#2a78d6"),
    ("thermal-const1.93e-18-x42", "thermal start, 10 mm in (x = 41.8 mm)", "#eb6834"),
    ("thermal-const1.93e-18", "thermal start, inlet channel exit (x = 31.8 mm)", "#1baf7a"),
    ("thermal-const1.93e-18-x66", "thermal start, 3/4 through (x = 65.9 mm)", "#eda100"),        # slot 4
    ("thermal-const1.93e-18-x75", "thermal start, 2 mm before outlet (x = 75.2 mm)", "#e87ba4"),  # slot 5
]
REF = ("edep", "500 eV injection, energy-dependent σ (previous model)", "#9a9994")   # neutral reference line
CONDS = CONDSETS["edep"]
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
M_BAF, EV = 156.325 * 1.66054e-27, 1.602177e-19
KT23 = 1.380649e-23 * 23 / EV
E_THERMAL = 0.1    # eV at box exit: below = thermalized (cold, gas-entrained ions sit at ~10-25 meV;
                   # ions that never cooled leave at ~1-400 eV)
COL_TH, COL_HOT = "#2a78d6", "#eb6834"   # reference palette slots 1-2


def split_hist(ax, x, hot, bins, n_injected):
    """Stacked histogram in % of ALL injected ions per bin (bars sum to the
    extraction fraction): thermalized ions (blue) under non-thermal (orange)."""
    xc = np.clip(x, bins[0], bins[-1])
    w = 100.0 / n_injected
    cnt, _, _ = ax.hist([xc[~hot], xc[hot]], bins=bins, stacked=True, color=[COL_TH, COL_HOT],
                        weights=[np.full((~hot).sum(), w), np.full(hot.sum(), w)],
                        edgecolor=SURF, lw=0.8, label=["thermalized", "non-thermal"])
    ax.set_ylim(0, max(w, np.max(cnt[-1])) * 1.45)   # headroom so the label never sits on a bar


def pct(k, n):
    p = 100.0 * k / n
    return "0%" if k == 0 else (f"{p:#.2g}%" if p < 10 else f"{p:.1f}%")


def split_legend(fig):
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(color=COL_TH, label=f"thermalized (< {E_THERMAL:g} eV at exit)"),
                        Patch(color=COL_HOT, label=f"non-thermal (≥ {E_THERMAL:g} eV)")],
               loc="upper right", ncol=2, frameon=False, fontsize=8, labelcolor=INK2)


def post_dir(f):
    return (f"{BASE}/fluor-cell-3d-ne-sealed/post" if f == 4 else f"{BASE}/fluor-scan/sccm{f}/post")


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan, np.nan
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return p, max(0.0, c - h), min(1.0, c + h)


def style(ax):
    ax.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--set", choices=sorted(CONDSETS), default="edep",
                    help="which cross-section conditions to plot (default edep: the original scan)")
    ap.add_argument("--out", default=None, help="output prefix (default results/figures/fluor3d_scan[_table])")
    a = ap.parse_args()
    global CONDS, FLOWS
    CONDS = CONDSETS[a.set]
    if a.set == "thermal":
        FLOWS = FLOWS_ALL
    primary = CONDS[0][0]
    ref = REF if a.set != "edep" else None
    if a.out is None:
        a.out = "results/figures/fluor3d_scan" + ("" if a.set == "edep" else f"_{a.set}")
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)

    rows, times, energies = [], {}, {}
    tags = [c[0] for c in CONDS] + ([ref[0]] if ref else [])
    for f in FLOWS:
        if ref and not glob.glob(f"{post_dir(f)}/tracer_{primary}*.out"):
            continue   # the table set only covers the flows it was run at
        for tag in tags:
            if ref and tag == ref[0] and f in EXCLUDE:
                continue
            files = sorted(glob.glob(f"{post_dir(f)}/tracer_{tag}.out") + glob.glob(f"{post_dir(f)}/tracer_{tag}_s*.out"))
            if not files:
                continue
            r = load(files)
            tr, lost, wall = classify(r)
            p, lo, hi = wilson(tr.sum(), len(r))
            rows.append(dict(sccm=f, cond=tag, n=len(r), transmitted=int(tr.sum()), lost=int(lost.sum()),
                             wall=int(wall.sum()), frac=p, lo=lo, hi=hi,
                             t_median_ms=float(np.median(r[tr, 11]) * 1e3) if tr.any() else np.nan,
                             E_median_eV=float(np.median(0.5 * M_BAF * (r[tr, 7:10] ** 2).sum(1) / EV)) if tr.any() else np.nan))
            if tag == primary:
                times[f] = r[tr, 11] * 1e3
                energies[f] = 0.5 * M_BAF * (r[tr, 7:10] ** 2).sum(1) / EV   # eV at box exit
    if not rows:
        sys.exit("no tracer output found")

    with open(a.out + "_summary.csv", "w") as fh:
        keys = list(rows[0])
        fh.write(",".join(keys) + "\n")
        for d in rows:
            fh.write(",".join(f"{d[k]:.6g}" if isinstance(d[k], float) else str(d[k]) for k in keys) + "\n")

    # --- extraction fraction vs flow
    fig, ax = plt.subplots(figsize=(7.5, 4.6), dpi=150)
    fig.patch.set_facecolor(SURF)
    style(ax)
    ends = []
    for tag, label, col in CONDS + ([ref] if ref else []):
        d = [x for x in rows if x["cond"] == tag]
        if not d:
            continue
        if ref and tag == ref[0]:   # previous model, drawn as a neutral dashed reference
            f = np.array([x["sccm"] for x in d]); p = np.array([x["frac"] for x in d]) * 100
            ax.plot(f, p, color=col, lw=1.5, ls="--", marker="o", ms=4, label=label, zorder=2)
            continue
        f = np.array([x["sccm"] for x in d]); p = np.array([x["frac"] for x in d]) * 100
        lo = np.array([x["lo"] for x in d]) * 100; hi = np.array([x["hi"] for x in d]) * 100
        ax.errorbar(f, p, yerr=[np.clip(p - lo, 0, None), np.clip(hi - p, 0, None)], color=col, lw=2,
                    marker="o", ms=6, capsize=3,
                    elinewidth=1, mec=SURF, mew=1.5, label=label, zorder=3)
        ends.append((p[-1], f[-1], label.split(" (")[0]))
        if tag == primary and a.set != "thermal":   # same condition, counting only ions that left thermalized
            k = np.array([(energies[x["sccm"]] < E_THERMAL).sum() for x in d])
            n = np.array([x["n"] for x in d])
            w = np.array([wilson(kk, nn) for kk, nn in zip(k, n)]) * 100
            ax.errorbar(f, w[:, 0], yerr=[w[:, 0] - w[:, 1], w[:, 2] - w[:, 0]], color=col, lw=1.5, ls="--",
                        marker="o", ms=5, mfc=SURF, mec=col, mew=1.2, capsize=2, elinewidth=0.8,
                        label=f"{label}, thermalized only (< {E_THERMAL:g} eV)", zorder=3)
    # direct labels at each series' last point, nudged apart where series end together
    ax.set_xscale("log")
    plotted = sorted({x["sccm"] for x in rows})
    xlo, xhi = (0.8, 170) if not ref else (plotted[0] / 1.25, plotted[-1] * 1.35)
    ax.set_xlim(xlo, xhi)
    fig.canvas.draw()
    placed = []
    px = fig.dpi / 72   # display pixels per point: positions are in pixels, offsets in points
    for y, x, text in sorted(ends):
        ypix = ax.transData.transform((x, y))[1]
        off = 0.0
        while any(abs(ypix + off - q) < 11 * px for q in placed):
            off += 11 * px
        placed.append(ypix + off)
        ax.annotate(text, (x, y), xytext=(6, off / px), textcoords="offset points", va="center", fontsize=8, color=INK2)
    ticks = ([t for t in (1, 2, 4, 8, 15, 20, 30, 58) if FLOWS[0] <= t <= FLOWS[-1]] if not ref else plotted)
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(t) for t in ticks])
    ax.minorticks_off()
    ax.set_xlabel("Buffer-gas (Ne) flow  [SCCM]", color=INK, fontsize=10)
    ax.set_ylabel("BaF⁺ extracted through outlet  [%]", color=INK, fontsize=10)
    ax.set_ylim(bottom=0)
    if UNDER_RESOLVED:   # one band: every flow >= 10 SCCM is under-resolved
        ax.axvspan(min(UNDER_RESOLVED) / 1.08, xhi, color=GRID, alpha=0.5, lw=0, zorder=0)
    ax.set_title("BaF⁺ extraction fraction vs. flow (" + ("thermal injection test" if a.set == "thermal" else "500 eV injection") + ", sealed cell)", color=INK, fontsize=11, loc="left",
                 pad=12 + 13 * int(np.ceil(len(ax.get_legend_handles_labels()[0]) / 2)))
    ax.legend(frameon=False, fontsize=8, labelcolor=INK2, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2,
              handlelength=1.8, columnspacing=1.2)
    fig.text(0.01, 0.01, "Error bars: 95% Wilson interval. Shaded (>= 10 SCCM): chamber grid coarser than twice the Ne\n"
             "mean free path, so port flow is approximate. " +
             ("1-2 SCCM extraction is hot ions flying straight through." if not ref else
              "Ions start as a 23 K Maxwellian (no beam) with constant\nσ = 1.93e-18 m², on the port axis at five depths."
              if a.set == "thermal" else
              "Dashed blue: ions leaving thermalized\n(< 0.1 eV); the rest at low flow are hot ions coasting straight through."),
             fontsize=7.5, color=INK2)
    fig.tight_layout(rect=(0, 0.08 if ref else 0.06, 1, 1))
    fig.savefig(a.out + "_extraction.png", facecolor=SURF)
    plt.close(fig)

    # --- pump-out time histograms, one panel per flow (thermalized / non-thermal stacked)
    fl = [f for f in FLOWS if f in times]
    ncol = min(4, len(fl)); nrow = int(np.ceil(len(fl) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(max(3.2 * ncol, 6.4), 2.6 * nrow + 0.4), dpi=150, squeeze=False)
    fig.patch.set_facecolor(SURF)
    for ax in axs.flat[len(fl):]:
        ax.set_visible(False)
    # times span us (ballistic transit at low flow) to tens of ms (diffusive
    # extraction at high flow), so one shared log-spaced axis for every panel
    bins = np.logspace(-3, 3, 49)   # 1 us .. 1 s, 8 bins/decade
    for ax, f in zip(axs.flat, fl):
        style(ax)
        t = times[f]
        n_all = next(x["n"] for x in rows if x["sccm"] == f and x["cond"] == primary)
        if len(t):
            hot = energies[f] >= E_THERMAL
            split_hist(ax, t, hot, bins, n_all)
            ax.axvline(np.median(t), color=INK, lw=1, ls="--")
            med = np.median(t)
            ax.text(0.97, 0.95, f"{pct(len(t), n_all)} extracted\n({pct((~hot).sum(), n_all)} thermalized)\nmedian " +
                    (f"{med*1e3:.1f} µs" if med < 1 else f"{med:.1f} ms"),
                    transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color=INK2,
                    bbox=dict(boxstyle="round,pad=0.25", fc=SURF, ec="none", alpha=0.9), zorder=5)
        else:
            ax.text(0.5, 0.5, f"0% extracted (of {n_all})", transform=ax.transAxes, ha="center", fontsize=8, color=INK2)
        ax.set_xscale("log")
        ax.set_xlim(bins[0], bins[-1])
        ax.set_xticks([1e-3, 1e-2, 1e-1, 1, 10, 100, 1000])
        ax.set_xticklabels(["1µs", "10µs", "0.1ms", "1ms", "10ms", "0.1s", "1s"])
        ax.minorticks_off()
        ax.set_title(f"{f} SCCM", fontsize=9, color=INK, loc="left")
        ax.set_xlabel("time from injection to exiting the outlet", fontsize=8, color=INK2)
        ax.set_ylabel("% of injected ions", fontsize=8, color=INK2)
    fig.suptitle(f"BaF⁺ pump-out time, injection to outlet ({CONDS[0][1]})", color=INK, fontsize=11, x=0.01, ha="left")
    split_legend(fig)
    fig.text(0.01, 0.01, "Dashed: median of all extracted ions.", fontsize=7.5, color=INK2)
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    fig.savefig(a.out + "_pumpout.png", facecolor=SURF)
    plt.close(fig)
    # --- kinetic energy at extraction, one panel per flow (thermalized / non-thermal stacked)
    fig, axs = plt.subplots(nrow, ncol, figsize=(max(3.2 * ncol, 6.4), 2.6 * nrow + 0.4), dpi=150, squeeze=False)
    fig.patch.set_facecolor(SURF)
    for ax in axs.flat[len(fl):]:
        ax.set_visible(False)
    ebins = np.logspace(-3, np.log10(600), 45)   # 1 meV .. 600 eV (injection is 500 eV)
    for ax, f in zip(axs.flat, fl):
        style(ax)
        E = energies[f]
        if len(E):
            hot = E >= E_THERMAL
            n_all = next(x["n"] for x in rows if x["sccm"] == f and x["cond"] == primary)
            split_hist(ax, E, hot, ebins, n_all)
            med = np.median(E)
            ax.axvline(med, color=INK, lw=1, ls="--")
            ax.axvline(1.5 * KT23, color=INK2, lw=1, ls=":")
            ax.text(0.97, 0.95, f"{pct(len(E), n_all)} extracted\n({pct((~hot).sum(), n_all)} thermalized)\nmedian " +
                    (f"{med*1e3:.1f} meV" if med < 1 else f"{med:.3g} eV"),
                    transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color=INK2,
                    bbox=dict(boxstyle="round,pad=0.25", fc=SURF, ec="none", alpha=0.9), zorder=5)
        else:
            ax.text(0.5, 0.5, "0% extracted", transform=ax.transAxes, ha="center", fontsize=8, color=INK2)
        ax.set_xscale("log")
        ax.set_xlim(ebins[0], ebins[-1])
        ax.set_xticks([1e-3, 1e-2, 1e-1, 1, 10, 100])
        ax.set_xticklabels(["1meV", "10meV", "0.1eV", "1eV", "10eV", "100eV"])
        ax.minorticks_off()
        ax.set_title(f"{f} SCCM", fontsize=9, color=INK, loc="left")
        ax.set_xlabel("kinetic energy at box exit", fontsize=8, color=INK2)
        ax.set_ylabel("% of injected ions", fontsize=8, color=INK2)
    fig.suptitle(f"Extracted BaF⁺ kinetic energy ({CONDS[0][1]}; injected at 500 eV)", color=INK, fontsize=11, x=0.01, ha="left")
    split_legend(fig)
    fig.text(0.01, 0.01, "Dashed: median of all extracted ions.  Dotted: 3/2 kT at 23 K (3.0 meV).", fontsize=7.5, color=INK2)
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    fig.savefig(a.out + "_energy.png", facecolor=SURF)
    plt.close(fig)
    print(f"wrote {a.out}_extraction.png, {a.out}_pumpout.png, {a.out}_energy.png, {a.out}_summary.csv")


if __name__ == "__main__":
    main()
