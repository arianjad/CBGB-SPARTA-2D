#!/usr/bin/env python3
"""Strict DSMC grid-rule check on a frozen (time-averaged) field of the plain-cylinder case.

For every cell with gas: mean free path and mean collision time from the VHS
parameters in ne.vss (the formulas of SPARTA's compute lambda/grid, single
species), the cell Knudsen number Kn = lambda / mean cell side, and the ratio
tau / DT. The rule (Steinmacher et al., arXiv:2609.11544) is Kn >= 1 and
DT < tau everywhere; here Kn >= KNMIN (default 0.5, cells up to 2 lambda) and DT < tau. Also prints the chamber-average density next to the
molecular-flow aperture law, n = 4 Q / (c_bar A), their eq. 2.
With --cellfrac (cellfrac.grid, written by in.ne_cyl3d_run after adaptation),
cells cut by the wall to a gas-volume fraction < --fmin are skipped, as in the
run's own check (their few-particle samples are not meaningful).

Usage: python tools/cyl3d_check_rule.py FROZEN_GRID --vss ne.vss --dt DT --sccm F [--cellfrac F --fmin 0.5]
"""
import argparse
import math

import numpy as np
import pandas as pd

KB = 1.380649e-23
SCCM = 4.4615e17                         # molecules/s per SCCM


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("frozen")
    ap.add_argument("--vss", required=True)
    ap.add_argument("--dt", type=float, required=True)
    ap.add_argument("--sccm", type=float, required=True)
    ap.add_argument("--cellfrac")
    ap.add_argument("--fmin", type=float, default=0.5)
    ap.add_argument("--aeff-mm2", type=float, default=2 * math.pi * 2.5**2,
                    help="port area for the aperture law, Clausing-weighted (default: two thin 5 mm apertures)")
    ap.add_argument("--knmin", type=float, default=0.5, help="cell <= lambda / KNMIN (0.5: up to 2 mean free paths)")
    ap.add_argument("--sid", default="Ne", help="species ID in the VSS file (Ne, He)")
    ap.add_argument("--mgas", type=float, default=20.1797, help="gas atomic mass (u)")
    ap.add_argument("--tgas", type=float, default=23.0, help="wall/inflow gas temperature (K), for the aperture law")
    ap.add_argument("--tfloor", type=float, default=3.0, help="floor on the sampled T (K), as TFLOOR in in.ne_cyl3d_run")
    a = ap.parse_args()
    M = a.mgas * 1.66053907e-27
    d, w, tref = [float(x) for x in next(l for l in open(a.vss) if l.split() and l.split()[0] == a.sid).split()[1:4]]
    v = pd.read_csv(a.frozen, sep=" ", header=None, skiprows=9, usecols=range(15), engine="c").values
    xc, yc, zc = v[:, 1], v[:, 2], v[:, 3]
    size = (v[:, 7:10] - v[:, 4:7]).mean(axis=1)
    n, T = v[:, 10], v[:, 14]
    gas = (n > 0) & (T > 0)
    if a.cellfrac:
        cf = pd.read_csv(a.cellfrac, sep=" ", header=None, skiprows=9, usecols=range(2), engine="c").values
        frac = dict(zip(cf[:, 0].astype(np.int64), cf[:, 1]))
        f = np.array([frac.get(int(i), 1.0) for i in v[:, 0]])
        print(f"skipping {(gas & (f < a.fmin)).sum()} wall-cut cells with gas fraction < {a.fmin}")
        gas &= f >= a.fmin
    T = np.maximum(T, a.tfloor)   # same floor as in.ne_cyl3d_run (TFLOOR): nearly empty cells sample T ~ 0
    lam = 1 / (math.pi * math.sqrt(2) * d**2 * n[gas] * (tref / T[gas])**(w - 0.5))
    tau = 1 / (2 * d**2 * n[gas] * math.sqrt(2 * math.pi * KB * tref / (M / 2)) * (T[gas] / tref)**(1 - w))
    kn = lam / size[gas]
    bad_kn, bad_dt = kn < a.knmin, tau <= a.dt   # knmin 0: fixed lattice, only dt < tau is enforced
    print(f"cells with gas {gas.sum()}; cell sizes {np.unique(np.round(size[gas]*1e3, 4))} mm")
    print(f"Kn_cell: min {kn.min():.3f}, 0.1% quantile {np.quantile(kn, 0.001):.3f}; cells with Kn < {a.knmin}: "
          f"{bad_kn.sum()} ({100*bad_kn.mean():.4f}%)")
    print(f"tau/DT: min {tau.min()/a.dt:.2f}; cells with DT >= tau: {bad_dt.sum()} ({100*bad_dt.mean():.4f}%)")
    if bad_kn.any():
        i = np.argsort(kn)[:5]
        for k in i:
            print(f"  worst Kn {kn[k]:.3f} at x,y,z = {xc[gas][k]*1e3:.2f}, {yc[gas][k]*1e3:.2f}, {zc[gas][k]*1e3:.2f} mm, "
                  f"n {n[gas][k]:.3e}, T {T[gas][k]:.1f} K, cell {size[gas][k]*1e3:.3f} mm")
    ch = gas & (xc > 0) & (xc < 0.04547) & (yc**2 + zc**2 < 0.01524**2)
    vol = np.prod(v[:, 7:10] - v[:, 4:7], axis=1)
    nch = (n * vol)[ch].sum() / vol[ch].sum()
    cbar = math.sqrt(8 * KB * a.tgas / (math.pi * M))
    neq2 = 4 * a.sccm * SCCM / (cbar * a.aeff_mm2 * 1e-6)
    print(f"chamber <n> {nch:.3e} m^-3 vs aperture law (eq. 2) {neq2:.3e} ({nch/neq2:.2f}x); "
          f"chamber <T> {(n*vol*T)[ch].sum()/(n*vol)[ch].sum():.1f} K")
    ok = not bad_kn.any() and not bad_dt.any()
    print("STRICT RULE " + ("SATISFIED" if ok else "VIOLATED"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
