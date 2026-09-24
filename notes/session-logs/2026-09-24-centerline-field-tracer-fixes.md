# Centerline field and tracer fixes (2026-09-24)

The B5 deck now appends SPARTA `vol` to new `field.grid` snapshots. The B5
plotter accepts both 11-column historical and 12-column new dumps. New dumps
show positive-flow-area zero-density cells in gray and include their zeros in
density profiles and region means. Historical zero-density cells remain
unclassified because those dumps did not record flow area.

The 2D tracer now rejects a positive-density cell with nonpositive saved
temperature by default. Explicit `--keep-unsampled` / `KEEP_UNSAMPLED=1`
borrows only temperature from the nearest positive-density cell with measured
temperature; density and velocity remain local. The standalone manifest records
the choice, and stderr records how many temperatures were borrowed. Combined
temperature plots exclude zero thermal sentinels from their denominator.

Checks passed in this worktree:

- `tools/test_field_io.py` including open-zero/solid geometry and zero-T windows.
- `tracer/test_partner_allocations.jl` with Julia 1.9.4, including production
  `build_field` interpolation for a positive-density zero-T centerline cell.
- `tools/test_standalone.sh`, including opt-in command and manifest recording.
- Full portable `tools/test.sh` suite, including geometry regeneration.
- `bash -n tools/run_2d_standalone.sh` and `git diff --check`.
- A one-step, zero-particle syntax smoke run of the revised full B5 deck with
  installed SPARTA 27Aug2026. Its `field.grid` header ended in `vol`.

The smoke run only validates deck parsing and dump layout. It does not validate
helium physics, the numerical effect of the policy on any production run, or
Eric's field, which was not available in this checkout.
