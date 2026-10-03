#!/usr/bin/env python3
"""Cross sections used by the cyl-3d runs vs the older models -> results/figures/cross_sections.png

Left:  BaF+-Ne, as the tracer uses it, vs BaF+ lab energy (gas at rest, 23 K):
       the old energy-dependent 1/v model with its floor (--sigma-model edep:
       sigma = max(2e-18 v_ref/v_coll, 1.5e-19) m^2), the computed table now used
       (tracer/data/bafplus_ne_sigma_tracer.tsv, central with the low-high band),
       and the Langevin capture cross section (x1.105, pure 1/v) for reference.
Right: Ne-Ne, through the viscosity it implies: Bird's 273 K VHS parameters
       extrapolated down, the ab initio values (quantum Chapman-Enskog on the
       Hellmann-Bich-Vogel 2008 potential, ~/orcd/scratch/ne_collision_integrals/),
       and the 5-60 K VSS fit now in cases/cyl-3d/ne.vss. Viscosity scales as
       1/(collision cross section), so lower eta = larger cross section.
Usage: python tools/plot_cross_sections.py [--out results/figures/cross_sections.png]
"""
import argparse
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

KB, AMU, EV = 1.380649e-23, 1.66053907e-27, 1.602176634e-19
M_NE, M_ION = 20.1797, 156.325
MU = M_NE * M_ION / (M_NE + M_ION)
C1, C2, C3, INK, MUTED = "#2a78d6", "#eb6834", "#1baf7a", "#1f1f1e", "#6b6a66"
# ab initio dilute-gas neon viscosity (uPa s), computed 2026-10-02 (see docstring)
ETA_AI = {5: 1.036, 7: 1.365, 10: 1.772, 15: 2.434, 20: 3.154, 23: 3.610, 25: 3.920, 30: 4.708, 40: 6.291,
          50: 7.825, 60: 9.285, 100: 14.40, 200: 24.12, 273.15: 29.90}


def vhs_eta(T, d, w, tref, alpha=1.0):
    m = M_NE * AMU
    eta_ref = 5 * (alpha + 1) * (alpha + 2) * math.sqrt(m * KB * tref / math.pi) / (
        4 * alpha * (5 - 2 * w) * (7 - 2 * w) * d**2)
    return eta_ref * (np.asarray(T) / tref) ** w


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/figures/cross_sections.png")
    a = ap.parse_args()
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))

    # BaF+-Ne vs lab energy
    Elab = np.logspace(-3, 4, 400)                      # eV
    vref = math.sqrt(8 * KB * 23.0 / (math.pi * M_NE * AMU))
    vrel = np.sqrt(2 * Elab * EV / (M_ION * AMU))
    vcoll = np.sqrt(vref**2 + vrel**2)
    edep = np.maximum(2e-18 * vref / vcoll, 1.5e-19)
    d = np.loadtxt("tracer/data/bafplus_ne_sigma_tracer.tsv", comments="#")
    Ecm_K = MU * AMU * vcoll**2 / 2 / KB
    tab = {k: np.exp(np.interp(np.log(Ecm_K), np.log(d[:, 0]), np.log(d[:, c]))) for k, c in (("c", 1), ("lo", 2), ("hi", 3))}
    alpha_ne = 0.3943e-30 * 4 * math.pi * 8.8541878e-12        # SI polarizability
    Ecm_J = MU * AMU * vrel**2 / 2
    langevin = 1.105 * math.pi * np.sqrt(2 * alpha_ne * EV**2 / (4 * math.pi * 8.8541878e-12)**2 / 2 / np.maximum(Ecm_J, 1e-30))
    ax[0].fill_between(Elab, tab["lo"], tab["hi"], color=C1, alpha=0.18, lw=0)
    ax[0].loglog(Elab, tab["c"], color=C1, lw=2.2, label="computed table (now used); band = low/high wall")
    ax[0].loglog(Elab, edep, color=C2, lw=2.2, label="old 1/v with 1.5e-19 floor (edep)")
    ax[0].loglog(Elab, langevin, color=MUTED, lw=1.2, ls="--", label="Langevin ×1.105 (pure 1/v)")
    ax[0].axvspan(250, 5000, color="#b9b8b2", alpha=0.25, lw=0)
    ax[0].text(330, 3e-21, "beam\nenergies", fontsize=8, color=MUTED)
    ax[0].set_xlabel("BaF⁺ lab energy (eV), Ne at 23 K"); ax[0].set_ylabel("collision cross section (m²)")
    ax[0].set_title("BaF⁺–Ne (ion tracer)", fontsize=11); ax[0].set_ylim(1e-21, 2e-17)
    ax[0].legend(fontsize=8, frameon=False, loc="lower left"); ax[0].grid(alpha=0.25, which="both")

    # Ne-Ne viscosity
    T = np.logspace(np.log10(4), np.log10(300), 200)
    ax[1].loglog(T, vhs_eta(T, 2.766e-10, 0.66, 273.15) * 1e6, color=C2, lw=2.2,
                 label="Bird VHS (d 2.77 Å, ω 0.66, 273 K): old runs")
    ax[1].loglog(T, vhs_eta(T, 4.59e-10, 0.883, 23.0, 1.43) * 1e6, color=C1, lw=2.2,
                 label="5–60 K VSS fit (d 4.59 Å, ω 0.883, α 1.43): cyl-3d runs")
    Ta = np.array(sorted(ETA_AI)); ax[1].loglog(Ta, [ETA_AI[t] for t in Ta], "o", color=C3, ms=5,
                                                label="ab initio (quantum Chapman–Enskog)")
    ax[1].axvline(23, color=MUTED, lw=1, ls=":"); ax[1].text(24, 1.1, "23 K cell", fontsize=8, color=MUTED)
    ax[1].annotate("+61% at 23 K", (23, vhs_eta(23, 2.766e-10, 0.66, 273.15) * 1e6), xytext=(8, 8),
                   textcoords="offset points", fontsize=8, color=INK)
    ax[1].set_xlabel("temperature (K)"); ax[1].set_ylabel("Ne viscosity (µPa·s)  [∝ 1/cross section]")
    ax[1].set_title("Ne–Ne (neon DSMC)", fontsize=11)
    ax[1].legend(fontsize=8, frameon=False, loc="upper left"); ax[1].grid(alpha=0.25, which="both")
    fig.tight_layout(); fig.savefig(a.out, dpi=150)
    print(f"-> {a.out}")


if __name__ == "__main__":
    main()
