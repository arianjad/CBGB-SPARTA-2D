# Plain cylinder with old-cell ports, 23 K neon: results (2026-10-06)

The cyl-3d `ports` scan uses:

- **Geometry:** a 30.5 mm bore × 45.5 mm chamber with the fluor-cell-3d ports
  (2.54 mm inlet through a 1.27 mm wall; 3.175 mm × 5.33 mm outlet channel). The rest
  of the old CAD features are removed.
- **Scan:** 5–50 SCCM neon at 23 K, BaF⁺ injected on axis at the inlet at
  0.25–5 keV in 250 eV steps, 2000 ions per point.
- **Cross sections:** the table-central BaF⁺–Ne cross section and the refit Ne–Ne VSS
  (`cases/cyl-3d/ne.vss`).
- **DSMC:** a fixed 0.4 mm lattice, with the time step below the local collision time
  everywhere. Every flow reached steady state (0.5% per 20 ms).

Dashboard: https://claude.ai/artifact/2Qavs8GntHojgLd1dwYzBY (Ne cylinder tab).
Figures and data: `results/figures/cyl3d_ports_*` (summary CSV: `cyl3d_ports_summary.csv`).

## Coverage: 190 of 200 points

| Missing | Why |
|---|---|
| 50 SCCM, 3500–5000 eV (7 points) | 3750–5000 eV jobs cancelled at the user's request (2026-10-05 19:45); the 3500 eV job hit its 24 h limit |
| 30 SCCM 250 eV; 40 SCCM 4250 and 4750 eV | Crashed with a Julia bus error on 2026-10-04 (shared compile-cache rewrite). Each kept its histograms and ion log, but only 1801 of 2000 final ion states (the last-finishing ions were lost, so the sample is biased). Their rescue jobs were cancelled on 2026-10-05 by mistake: the audit counted `_xr.csv` files as finished points instead of `.out` files |

The two 50 SCCM GIFs at 3000 and 3250 eV were never rendered (time limit). The
partial `.out.tmp` files are kept.

## Main results

- **Best point:** thermalized (extracted, < 0.1 eV) is **36.8%** (95% interval
  34.7–38.9) at 40 SCCM, 5 keV.
- **Next best:** 34.5% (35 SCCM, 4.75 keV), 34.4% (40 SCCM, 4.5 keV) and 33.1%
  (30 SCCM, 5 keV).
- **Extracted vs thermalized:** above ~15 SCCM essentially every extracted ion is
  thermalized. The ~0.03 eV exit energy is the gas flow carrying the ions, not heat (see
  `docs/sanity-checks-cyl3d.md`).
- **Trends:** the best energy per flow rises with flow and is at or near the top of the
  range from 30 SCCM up:

  | Flow (SCCM) | Best energy | Thermalized |
  |---|---|---|
  | 20 | 3.75 keV | 21.9% |
  | 25 | 4.5 keV | 30.6% |
  | 30 | 5 keV | 33.1% |
  | 35 | 4.75 keV | 34.5% |
  | 40 | 5 keV | 36.8% |
  | 45 | 4.75 keV | 32.8% |

  The optimum may lie above 5 keV at 30–45 SCCM. 50 SCCM peaks at 18.9% (3.25 keV),
  but its 3.5–5 keV points are missing.
- **Low energies at high flow lose ions out the inlet:** at 500 eV, the fraction
  leaving back through the inlet rises with flow:

  | Flow (SCCM) | Back out the inlet |
  |---|---|
  | 20 | 11% |
  | 30 | 50% |
  | 40 | 82% |
  | 50 | 95% |

  Ions stop close to the inlet, where gas streams back out through it.

## Comparison with the old fluor-cell-3d (table-central)

The old cell (CAD geometry, the same ports) gave a chamber density of
1.25e21 m⁻³ per SCCM; this cylinder gives 1.18–1.28e21. So the same flow is a matched
density within 5–8% (the cylinder is slightly thinner). In the table, "old" is the old
fluor-cell scan (`results/figures/fluor3d_scan_kev_summary.csv`; 1500–12000 ions per
point) and "new" is this scan. Values are extracted %, with thermalized % after the
slash (the old scan's thermalized at 500 eV equals extracted).

| SCCM | n old / new (m⁻³) | 0.5 keV, old | 0.5 keV, new | 2 keV, old | 2 keV, new | 5 keV, old | 5 keV, new |
|---|---|---|---|---|---|---|---|
| 5 | 6.9e21 / 6.4e21 | 3.2 | 3.6 / 0.9 | 12.5 / 0.0 | 14.4 / 0.0 | 28.6 / 0.0 | 35.9 / 0.0 |
| 10 | 1.30e22 / 1.21e22 | 4.6 | 5.3 | 7.8 / 5.9 | 7.3 / 5.3 | 10.5 / 1.5 | 11.2 / 1.6 |
| 15 | 1.90e22 / 1.77e22 | 6.5 | 8.1 | 12.9 / 12.7 | 14.8 | 12.5 / 10.4 | 12.4 / 9.6 |
| 20 | 2.52e22 / 2.34e22 | 8.6 | 8.1 | 15.2 | 17.1 | 20.2 / 19.4 | 21.1 / 20.2 |
| 25 | 3.13e22 / 2.92e22 | 8.6 | 8.1 | 14.5 | 21.2 | 22.0 | 29.6 |
| 30 | 3.74e22 / 3.50e22 | 5.6 | 4.7 | 14.7 | 17.9 | 22.4 | 33.1 |
| 40 | 4.98e22 / 4.69e22 | – | 1.9 | 13.2 | 17.5 | 25.8 | 36.8 |
| 45 | 5.65e22 / 5.27e22 | – | 0.9 | 13.2 | 13.2 | 23.0 | 32.2 |
| 50 | 6.24e22 / 5.88e22 | – | 0.4 | 10.8 | 10.2 | 22.6 | – |

**Reading the comparison:**

- **Up to 20 SCCM the two cells agree within ~2 points** at every energy. Statistical
  uncertainty is ±0.6–1.5 points, so the plain cylinder reproduces the old cell where
  ions thermalize mid-chamber.
- **From 25 SCCM up, at 2–5 keV, the plain cylinder extracts more:** +7 to +11 points
  at 5 keV (for example 33% vs 22% at 30 SCCM, 37% vs 26% at 40 SCCM), and +3 to
  +7 points at 2 keV. At 0.5 keV, and at 45–50 SCCM and 2 keV, they agree again.
- **The two runs differ in more than geometry,** so the high-flow gap can't be pinned
  on the CAD features alone:
  1. **Ne–Ne model:** the old cell used Bird's 273 K VHS parameters
     (`Ne 2.766e-10 0.66 273.15`), 61% too viscous at 23 K. The cylinder uses the
     ab initio refit. Viscosity sets the flow pattern that carries thermal ions to the
     outlet.
  2. **Geometry:** internal features of the old cell versus a smooth bore.
  3. **Density:** at the same flow, the cylinder is ~7% thinner.
  4. **Lattice:** the cylinder's lattice is a fixed 0.4 mm, about 2.4–22 lattice
     spacings per mean free path (see `docs/TODO.md`); the old grids differed.

  Item 1 can be separated by rerunning one high-flow cylinder field with Bird's
  parameters.

## Known biases common to both scans

- **Start-of-flight free-path rule.** Each free path was drawn with the gas density at
  the flight's start point. Ions starting in the low-density inlet bore go 2–4 mm
  (12–22%) too deep, and the entrance dE/dx is ~25% low. Extraction is probably
  somewhat overestimated, by the same rule in both scans.
  - **Fix:** `--max-flight`, on for new runs, off here for consistency
    (`docs/sanity-checks-cyl3d.md`).
- **BaF⁺ treated as Ba⁺ for collisions:** the F atom is ignored (about +16–22% stopping
  at keV would come from it), and collision-induced dissociation above ~56 eV lab is
  not modelled.
- **Coarse fixed lattice at high flow:** the lattice-convergence test is still on the
  TODO list.
