# Saved-field reader migration correction

Implemented and independently reviewed on 11 November 2026. Baseline: `077f30ef8e4b0f2a33a1d1cbf67f707023d53e77`. Historical reference: DSMC-sims `fix/centerline-analysis` at `06669024e51edd662f176428e1be1c52b011499f`; only evidenced reader behavior was ported, with no runtime import from another repository.

- Preserve time-mean density and density-weighted saved-window velocity.
- Exclude zero/NaN temperature sentinels from measured temperature weights, keeping zero when no measured temperature exists.
- Accept consistent 11-column legacy frames or 12-column frames with trailing `vol` flow area.
- Include zero-density open cells in area-weighted regional density summaries.
- Validate IDs, finite/stable geometry, and finite/nonnegative/stable area across every eligible trend frame, not just the averaging tail.

The two-window regression changes temperature from 7.5 to 10 while preserving velocity 17.5. Six default-fraction grid/area cases failed before review correction and now reject invalid input. All 16 B5 regressions and existing complete-frame IO tests pass, independently rerun. No solver runs, scientific defaults, dependencies, historical outputs or unrelated workflows changed.

Run `tools/test_field_reader_regressions.py` and `tools/test_field_io.py` with the existing analysis Python and BLAS/OMP threads set to 1. On this Windows installation, that environment's `Library/bin` is needed in the process-local PATH for native plotting. No persistent setting change is required.

Detailed baseline/final evidence and the saved-data rendering check are retained in the companion SPARTA project at `reports/2026-11-11/reader-corrections/REPORT.md`; this is documentation provenance, not a runtime dependency. These checks cover saved-data interpretation, not gas equilibrium, solver accuracy or convergence. Saved-window temperature is a density-weighted summary, not pooled variance.

Glossary: a sentinel denotes an unmeasured value; flow area is the portion of a cell available to gas; a regression is a small test that fails before a correction and passes after it.
