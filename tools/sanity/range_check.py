#!/usr/bin/env python3
"""Hand-calculation sanity check of the cyl-3d ports tracer runs: stopping depth,
collisions to thermalize and the dE/dx curve, predicted from the cross-section
table alone, vs what the tracer produced.

Physics (same collision model as tracer/ParticleTracing3D.jl, nothing else taken
from the tracer): BaF+ (m1 = 156.325 u) hits Ne (m2 = 20.18 u) with rate
n sigma(E_cm) v_coll and scatters isotropically in the centre-of-mass frame.
For a target at rest, one collision leaves the ion with the fraction
f = 1 - g u of its energy, g = 4 m1 m2/(m1+m2)^2 = 0.405, u uniform on [0, 1]:
  mean energy loss       <dE> = (g/2) E = 0.2025 E  -> dE/ds = n sigma(E) (g/2) E
  typical (log-mean) ion <ln f> = -0.2373          -> N(E0 -> 1 eV) = ln(E0/1 eV)/0.2373
  forward velocity       <v_x'> = v_x m1/(m1+m2)   -> direction cosine shrinks by
                         (m1/M) / exp(<ln f>/2) = 0.997 per collision (projected vs path length)
Free path between collisions 1/(n sigma). The neon density n(x) is read along the
beam axis from the frozen DSMC field (cells with r < 0.6 mm); the run's chamber
average is shown for comparison. Thermal motion of the 23 K neon is kept in
v_coll (it matters only below ~1 eV).

Sim side: <tag>_ions.csv (x where each ion first drops below 1 eV, collisions so
far) and <tag>_xr.csv (energy deposited per x bin, summed over r).

Reads the on-axis cells of each flow's frozen.grid; set AXIS_DIR to a folder of axis_<F>.txt files
(awk 'NR>9 && ($3*$3+$4*$4)<3.6e-7 {print $2,$11,$15}' frozen.grid) to skip parsing the full grid.
Writes results/sanity/range_check.txt and results/sanity/range_check.png.
"""
import os
import sys

import numpy as np
import pandas as pd

BASE = os.path.expanduser("~/orcd/scratch/cbgb-sparta-3d/cyl-ports")
TABLE = "tracer/data/bafplus_ne_sigma_tracer.tsv"
CASES = [(10, 500), (15, 2000), (20, 1000), (30, 5000), (40, 3000)]   # (SCCM, eV)
X0 = -1.27e-3          # beam start: outer face of the inlet aperture (m)
OUT = "results/sanity"

kB, amu, eV = 1.380649e-23, 1.66053906660e-27, 1.602176634e-19
m1, m2 = 156.325, 20.1797
M = m1 + m2
G = 4 * m1 * m2 / M**2                       # max fractional energy loss per collision
u = (np.arange(200000) + 0.5) / 200000
LNF = np.log(1 - G * u)
LNF_MEAN, LNF_SD = LNF.mean(), LNF.std()     # per-collision log energy fraction
COS_DECAY = (m1 / M) / np.exp(LNF_MEAN / 2)  # forward-velocity loss over speed loss (< 1)
MU = m1 * m2 / M * amu

tab = np.array([l.split() for l in open(TABLE) if l.strip() and not l.startswith("#")], float)
LE, LQ = np.log(tab[:, 0]), np.log(tab[:, 1])


def sigma(E_lab_eV, T):
    """Table sigma at the tracer's E_cm = mu v_coll^2/2, v_coll^2 = 8kT/(pi m2) + V^2."""
    V2 = 2 * E_lab_eV * eV / (m1 * amu)
    vc2 = 8 * kB * T / (np.pi * m2 * amu) + V2
    return np.exp(np.interp(np.log(MU * vc2 / 2 / kB), LE, LQ))


def axis_profile(f):
    """On-axis neon density and temperature: cells with r < 0.6 mm, averaged per x.
    Uses AXIS_DIR/axis_<f>.txt (columns x n T, cut from frozen.grid with awk) when given."""
    d = os.environ.get("AXIS_DIR")
    if d:
        g = pd.read_csv(f"{d}/axis_{f}.txt", sep=r"\s+", header=None, names=["x", "n", "T"])
    else:
        g = pd.read_csv(f"{BASE}/sccm{f}/post/frozen.grid", sep=r"\s+", skiprows=9, header=None,
                        usecols=[1, 2, 3, 10, 14], names=["x", "y", "z", "n", "T"])
        g = g[np.hypot(g.y, g.z) < 0.6e-3]
    p = g.groupby("x")[["n", "T"]].mean().reset_index()
    return p.x.values, p.n.values, p["T"].values      # (p.T would be the transpose)


def chamber_n(f):
    s = open(f"{BASE}/sccm{f}/post/chamber.txt").read()
    return float(s.split("<n>")[1].split()[0])


def march(E0, xs, ns, Ts, ds=2e-6, Emin=1.0):
    """Mean trajectory: returns arrays x, E_mean, dEdx_projected (eV/m), and the x where
    the typical (log-mean) ion reaches Emin, its collision count, and path length."""
    x, s, k, cos = X0, 0.0, 0.0, 1.0
    lnE = np.log(E0)                 # log-mean energy (typical ion)
    Em = E0                          # mean energy
    out_x, out_E, out_d = [], [], []
    x_stop = None
    while x < 0.050 and Em > 1e-3:
        n = np.interp(x, xs, ns, left=0.0, right=0.0)
        T = np.interp(x, xs, Ts)
        dk = n * sigma(np.exp(lnE), T) * ds
        dEds = n * sigma(Em, T) * (G / 2) * Em            # mean energy loss per path
        out_x.append(x); out_E.append(Em); out_d.append(dEds / max(cos, 1e-3))
        Em -= dEds * ds
        lnE += LNF_MEAN * dk
        k += dk; s += ds
        x += cos * ds
        cos *= COS_DECAY ** dk
        if x_stop is None and lnE <= np.log(Emin):
            x_stop, k_stop, s_stop, cos_stop = x, k, s, cos
        if x_stop is not None and Em < Emin:
            break
    if x_stop is None:
        x_stop = k_stop = s_stop = cos_stop = np.nan
    return np.array(out_x), np.array(out_E), np.array(out_d), x_stop, k_stop, s_stop, cos_stop


def sim_stats(f, E):
    tag = f"{BASE}/sccm{f}/post/tracer_table-central-e{E}_s1"
    ions = pd.read_csv(f"{tag}_ions.csv", comment="#")
    ok = ions.x_1eV.notna()
    xr = pd.read_csv(f"{tag}_xr.csv", comment="#")
    hdr = open(f"{tag}_xr.csv").readline()
    dx = float(hdr.split("dx=")[1].split()[0])
    xlo = float(hdr.split("xlo=")[1].split()[0])
    N = len(ions)
    ed = xr.groupby("ix").edep_eV.sum()
    xc = xlo + (ed.index.values - 0.5) * dx
    dedx = ed.values / N / dx                                # eV per m per ion
    return dict(N=N, frac=ok.mean(), x_mean=ions.x_1eV[ok].mean(), x_med=ions.x_1eV[ok].median(),
                x_sd=ions.x_1eV[ok].std(), k_mean=ions.n_1eV[ok].mean(), k_sd=ions.n_1eV[ok].std(),
                xc=xc, dedx=dedx, x_peak=xc[np.argmax(dedx)], d_peak=dedx.max())


def main():
    os.makedirs(OUT, exist_ok=True)
    rows, curves = [], []
    for f, E in CASES:
        xs, ns, Ts = axis_profile(f)
        nch = chamber_n(f)
        px, pE, pd_, xs1, ks1, ss1, cs1 = march(E, xs, ns, Ts)
        # uniform chamber-average density, same model, for comparison
        ux, uE, ud, xu1, ku1, su1, cu1 = march(E, np.array([X0, 0.06]), np.array([nch, nch]),
                                             np.array([23.0, 23.0]))
        k_pred = np.log(E / 1.0) / -LNF_MEAN
        k_sd = np.sqrt(k_pred) * LNF_SD / -LNF_MEAN
        s = sim_stats(f, E)
        n_ax = np.interp(np.clip(px, 0, None), xs, ns).mean()
        rows.append((f, E, nch, n_ax, k_pred, k_sd, ks1, xs1, ss1, cs1, xu1, px[np.argmax(pd_)], pd_.max(), s))
        curves.append((f, E, px, pd_, ux, ud, s, xs1))
    with open(f"{OUT}/range_check.txt", "w") as fo:
        w = lambda t: (print(t), fo.write(t + "\n"))
        w("Hand calculation vs tracer, cyl-3d ports, table-central. x in mm, 0 = chamber side of the inlet wall;")
        w("ions start at x = -1.27 mm. 'calc' = mean-trajectory march through the on-axis DSMC density;")
        w("'uniform' = same with the chamber-average density. Sim: 2000 ions, mean over ions that reach 1 eV.")
        w(f"per collision: <dE>/E = {G/2:.4f}, <ln f> = {LNF_MEAN:.4f} (sd {LNF_SD:.4f}), direction cosine x{COS_DECAY:.4f}")
        w("")
        w(f"{'SCCM':>4} {'E/eV':>5} {'n_ch':>8} {'n_axis':>8} | {'N_1eV calc':>11} {'sim':>11} | "
          f"{'x_1eV calc':>10} {'uniform':>7} {'sim mean':>8} {'sim med':>7} {'sim sd':>6} {'reach%':>6} | "
          f"{'path calc':>9} {'cos':>4} | {'dE/dx peak x calc':>17} {'sim':>5} | {'peak eV/mm calc':>15} {'sim':>6}")
        for (f, E, nch, nax, kp, ksd, ks1, xs1, ss1, cs1, xu1, xpk, dpk, s) in rows:
            w(f"{f:4d} {E:5d} {nch:8.2e} {nax:8.2e} | {kp:5.1f}±{ksd:4.1f} {s['k_mean']:5.1f}±{s['k_sd']:4.1f} | "
              f"{xs1*1e3:10.1f} {xu1*1e3:7.1f} {s['x_mean']*1e3:8.1f} {s['x_med']*1e3:7.1f} {s['x_sd']*1e3:6.1f} "
              f"{100*s['frac']:6.1f} | {ss1*1e3:9.1f} {cs1:4.2f} | {xpk*1e3:17.2f} {s['x_peak']*1e3:5.2f} | "
              f"{dpk/1e3:15.1f} {s['d_peak']/1e3:6.1f}")
    plot(curves)


def plot(curves):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    SIM, CALC, UNI, INK, MUTED = "#2a78d6", "#eb6834", "#898781", "#1f1f1e", "#6b6a64"
    fig, axs = plt.subplots(1, len(curves), figsize=(3.3 * len(curves), 3.4), sharey=False)
    for ax, (f, E, px, pd_, ux, ud, s, xs1) in zip(axs, curves):
        ax.plot(s["xc"] * 1e3, s["dedx"] / 1e3, color=SIM, lw=2, label="tracer (2000 ions)")
        ax.plot(px * 1e3, pd_ / 1e3, color=CALC, lw=2, ls="--", label="hand calc, on-axis n(x)")
        ax.plot(ux * 1e3, ud / 1e3, color=UNI, lw=1.5, ls=":", label="hand calc, chamber-average n")
        ax.axvline(s["x_mean"] * 1e3, color=SIM, lw=1, alpha=0.6)
        ax.axvline(xs1 * 1e3, color=CALC, lw=1, ls="--", alpha=0.6)
        ax.set_title(f"{f} SCCM, {E/1000:g} keV", fontsize=10, color=INK)
        ax.set_xlim(-2, max(xs1 * 1e3 * 1.8, s["x_mean"] * 1e3 * 1.8))
        ax.set_xlabel("x (mm)", fontsize=9, color=MUTED)
        ax.tick_params(labelsize=8, colors=MUTED)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.grid(alpha=0.25, lw=0.5)
    axs[0].set_ylabel("energy loss per ion (eV/mm)", fontsize=9, color=MUTED)
    axs[0].legend(fontsize=7.5, frameon=False, loc="upper right")
    fig.suptitle("dE/dx vs depth: tracer vs hand calculation (vertical lines: mean depth where ions first reach 1 eV)",
                 fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(f"{OUT}/range_check.png", dpi=110)


if __name__ == "__main__":
    sys.exit(main())
