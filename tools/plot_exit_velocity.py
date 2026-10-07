#!/usr/bin/env python3
"""Exit-velocity analysis of the BaF+ ions that leave through the INLET port (cyl-3d).

For one tracer output (TAG.out + TAG_ions.csv from sbatch_cyl3d_tracer.sh):
  - transverse velocity components v_y, v_z of the ions that leave back through the inlet,
    as a histogram against a 23 K BaF+ Maxwellian (sigma = sqrt(kT/m)) and a Gaussian with
    the fitted transverse temperature T_perp = m <v_y^2 + v_z^2> / (2 k);
  - axial velocity v_x: mean drift and the spread about it (T_par = m var(v_x) / k);
  - how thermalized they are: fraction below 0.1 eV, and collisions after the ion first
    dropped below 0.1 eV (total collisions minus n_0p1eV from the ion log).
Only translation is modelled; rotational / vibrational states of BaF+ are not.
Ions are recorded where they leave the simulation box (x = -22.4 mm, in the plume
outside the inlet), so the velocities include whatever the expanding gas did to them there.

Usage: python tools/plot_exit_velocity.py POSTDIR/tracer_table-central-e250_s1 --out results/figures/NAME [--title ...]
Writes NAME.png and NAME.txt.
"""
import argparse
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tracer_classify3d import CYL_BOX, CYL_PORTS_OUTLET_X, classify, load  # noqa: E402

KB, AMU, EV = 1.380649e-23, 1.66053906660e-27, 1.602176634e-19
M = 156.325 * AMU
BLUE, ORANGE, INK, MUTED = "#2a78d6", "#eb6834", "#1f1f1e", "#6b6a66"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tag", help="path prefix of the tracer output (without .out)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default="")
    ap.add_argument("--tgas", type=float, default=23.0)
    a = ap.parse_args()
    r = load([a.tag + ".out"])
    tr, back, wall = classify(r, CYL_PORTS_OUTLET_X, CYL_BOX)
    v = r[:, 7:10]
    ke = 0.5 * M * (v**2).sum(1) / EV
    ions = pd.read_csv(a.tag + "_ions.csv", comment="#").set_index("idx")
    idx = r[:, 0].astype(int)
    n01 = ions.reindex(idx)["n_0p1eV"].to_numpy()
    after = r[:, 10] - n01                       # collisions after first dropping below 0.1 eV
    b = back
    vy, vz, vx = v[b, 1], v[b, 2], v[b, 0]
    Tperp = M * np.mean(vy**2 + vz**2) / (2 * KB)
    Tpar = M * np.var(vx) / KB
    # plume divergence: in an expanding jet the radial velocity grows with the distance from
    # the axis (v_r ~ k r). Remove that correlated part to get the random (thermal) part.
    yy, zz = r[b, 5], r[b, 6]
    rr = np.hypot(yy, zz) + 1e-12
    vr = (yy * vy + zz * vz) / rr
    vth = (yy * vz - zz * vy) / rr                       # azimuthal component: no divergence term
    kfit = np.sum(vr * rr) / np.sum(rr * rr)
    vr_res = vr - kfit * rr
    Trand = M * np.mean(vr_res**2 + vth**2) / (2 * KB)
    Tazi = M * np.mean(vth**2) / KB
    sig = math.sqrt(KB * a.tgas / M)
    lines = [
        f"{a.title}",
        f"ions traced {len(r)}: back out the inlet {b.sum()} ({100*b.mean():.1f}%), wall {wall.sum()} ({100*wall.mean():.1f}%), "
        f"other exit {tr.sum()}",
        f"back-out ions below 0.1 eV: {100*np.mean(ke[b] < 0.1):.1f}%; median KE {np.median(ke[b]):.4f} eV",
        f"collisions after first < 0.1 eV (back-out ions): median {np.nanmedian(after[b]):.0f}, "
        f"5th percentile {np.nanpercentile(after[b], 5):.0f}, minimum {np.nanmin(after[b]):.0f}",
        f"transverse temperature T_perp = {Tperp:.1f} K (gas {a.tgas:g} K); rms v_y {np.std(vy):.1f}, v_z {np.std(vz):.1f} m/s "
        f"(23 K Maxwellian: {sig:.1f} m/s)",
        f"plume divergence: v_r = {kfit:.0f} m/s per m x r (median exit radius {np.median(rr)*1e3:.1f} mm); "
        f"after removing it the transverse temperature is {Trand:.1f} K (azimuthal component alone: {Tazi:.1f} K)",
        f"axial: mean v_x {np.mean(vx):.0f} m/s (negative = back through the inlet), T_par about the drift = {Tpar:.1f} K",
    ]
    open(a.out + ".txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.2))
    lim = max(4 * sig, np.percentile(np.abs(np.r_[vy, vz]), 99.5))
    bins = np.linspace(-lim, lim, 41)
    ax[0].hist(np.r_[vy, vz], bins=bins, density=True, color=BLUE, alpha=0.55, label=f"v_y and v_z of back-out ions (n={b.sum()})")
    x = np.linspace(-lim, lim, 400)
    g = lambda s: np.exp(-x**2 / (2 * s**2)) / (s * math.sqrt(2 * math.pi))
    ax[0].plot(x, g(sig), color=INK, lw=2, ls="--", label=f"{a.tgas:g} K Maxwellian")
    ax[0].plot(x, g(math.sqrt(KB * Tperp / M)), color=ORANGE, lw=2, label=f"fit: T⊥ = {Tperp:.1f} K")
    ax[0].set_xlabel("transverse velocity component (m/s)"); ax[0].set_ylabel("probability density (s/m)")
    ax[0].set_title("Transverse velocity at exit", fontsize=11, color=INK)
    ax[0].legend(fontsize=8, frameon=False)
    ax[1].hist(vx, bins=40, color=BLUE, alpha=0.55)
    ax[1].axvline(np.mean(vx), color=ORANGE, lw=2, label=f"mean {np.mean(vx):.0f} m/s, T∥ = {Tpar:.1f} K about it")
    ax[1].set_xlabel("axial velocity v_x (m/s; negative = out through the inlet)"); ax[1].set_ylabel("ions")
    ax[1].set_title("Axial velocity at exit", fontsize=11, color=INK)
    ax[1].legend(fontsize=8, frameon=False)
    ax[2].scatter(rr * 1e3, vr, s=4, color=BLUE, alpha=0.4, lw=0, label="back-out ions")
    xx = np.linspace(0, rr.max(), 50)
    ax[2].plot(xx * 1e3, kfit * xx, color=ORANGE, lw=2, label=f"plume divergence v_r = {kfit*1e-3:.1f} (m/s)/mm × r")
    ax[2].set_xlabel("exit radius r (mm)"); ax[2].set_ylabel("radial velocity v_r (m/s)")
    ax[2].set_title(f"Divergence removed: random T⊥ = {Trand:.1f} K", fontsize=11, color=INK)
    ax[2].legend(fontsize=8, frameon=False)
    for x_ in ax:
        x_.grid(alpha=0.25)
        for sp in ("top", "right"):
            x_.spines[sp].set_visible(False)
    fig.suptitle(a.title or os.path.basename(a.tag), fontsize=11, color=INK)
    fig.tight_layout(); fig.savefig(a.out + ".png", dpi=130)
    print("->", a.out + ".png")


if __name__ == "__main__":
    main()
