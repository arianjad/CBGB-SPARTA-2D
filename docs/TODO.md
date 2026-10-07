# To do

## Open

- Tracer free-flight bias (found 2026-10-05, tools/sanity/). The tracer drew each
  free path with the gas density at the flight's start point, so ions starting in the
  low-density inlet bore went 2-4 mm (12-22%) too deep and the entrance dE/dx was ~25%
  low. Fix: `--max-flight` (null steps), wired into tools/wsl/sbatch_cyl3d_tracer.sh
  as MAXFLIGHT (0.2 mm for every run except RUNBASE cyl-ports, which keeps 0 so the
  current 200-point campaign is self-consistent). When that campaign is done: make
  0.2 mm the tracer's own default, and note the bias next to the cyl-ports results
  (and the old fluor-cell results, which used the same logic).
- Lattice-size convergence vs density (cyl-3d). Find which DSMC lattice
  spacings are acceptable at which neon density, i.e. the largest
  lattice / mean-free-path ratio at which BaF+ extraction (and stopping depth)
  stop changing.
  The ports scan uses a fixed 0.4 mm lattice. Lattice / mean free path
  goes from ~2.4 (5 SCCM, n = 6.4e21 m^-3) to ~22 (50 SCCM, n = 5.9e22 m^-3).
  The only test so far (old fluor-cell, 17 SCCM) found extraction unchanged
  between ~8 and ~4 mean free paths; nothing above ~8 has been checked, so the
  30-50 SCCM results are the least certain.
  Plan: at 2-3 densities (e.g. 20, 35, 50 SCCM), rerun the neon field at
  0.4, 0.2 and (where affordable) 0.1 mm and retrace a few energies (e.g. 1,
  3, 5 keV) with the same seeds; compare extraction, thermalized (extracted)
  and chamber density. Derive a rule of thumb (max lattice / mfp) and apply it
  to future runs.
- BaF+ F-atom stopping correction (optional, ~+20% stopping at keV):
  Bragg-additive F-Ne nuclear stopping on top of the Ba+-Ne table.
- Old fluor-cell items (from docs/HANDOFF-subMIT.md section 5): FINE-grid
  write-up at 17 SCCM, the 0.10 mm port grid, and the agreed sanity checks
  (diffusion vs Chapman-Enskog / Ba+ mobility, flow-following, range vs
  SRIM/LSS, chunk binomial variance, neon mass balance).

## Done

(move items here with the date and the commit that closed them)
