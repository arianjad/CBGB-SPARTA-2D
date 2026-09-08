# Validation

## Release check

Checked on 2026-09-08 in Ubuntu 22.04 under WSL2, with Python 3.10 and
Julia 1.9.4:

- A fresh build of the pinned SPARTA source passed all 242 upstream tests.
- The portable checks passed in a separate clean clone.
- The analytic-prefill helium example completed 120,000 steps (12 ms) on
  eight MPI ranks. Its grid had no face-adjacent size jumps above 2:1.
- A 1,000-molecule BaF run on four Julia threads used the final field at
  timestep 120000, retained 1,000 birth/crossing records and 20 complete
  recorded paths, and produced all five figures. Both run manifests were
  complete with zero return codes; exit and crossing scorers agreed on
  the escaped count. Helium and trajectory figures were visually inspected.

The README starts with four helium ranks for a smaller machine; set `RANKS`
to match available CPU capacity. These checks establish the workflow at the
settings above, not a universal runtime or precision guarantee.

A fresh Terra session independently followed the README from public snapshot
`6e9bb6c`, using a new Python environment, Julia installation/depot, and SPARTA
build. The four-rank 12 ms helium run and 1,000-molecule run both completed;
the latter produced 20 paths and five figures. Its reported sampler warnings
come from initialization of a table that exact-mode collisions do not use;
see [Gotchas](2d-gotchas.md#molecule-source).

The README trajectory preview was regenerated from the retained tutorial paths
as two signed Cartesian projections, without rerunning the simulation or
changing its recorded results. The overview and exit zoom were visually
inspected. Regression checks cover signed centerline crossings, both transverse
coordinates, terminal clipping fractions, wall fates, zero-collision paths,
and extracted-only figures. The portable suite passed after this plotting update.

## Portable checks

From the repository root after installation:

```bash
source .venv/bin/activate
bash tools/test.sh
```

These checks use synthetic data and temporary directories. They cover
clipped wall/domain contacts, exit statistics, first-crossing accounting,
trajectory coloring, geometry regeneration, and launcher refusal/failure
receipts. They also exercise overlapping independent launches, complete-frame
selection from growing dumps, frozen plot inputs, and mid-run molecule
launching up to the Julia boundary. They do not download packages or launch
SPARTA.

Additional installed-tracer checks on the release helium field passed for
`N=1` (1,203 collisions) and a three-particle, negligible-cross-section case
(zero collisions for every particle). Both passed actual leg accounting and
endpoint continuity. The shuffled-record diagnostic is no longer a runtime
requirement. These are software checks, not physical estimates from such
small samples.

The updated standalone wrapper also completed an `N=1`, `TRAJPRINT=1`,
`UNTIL_MS=5` run using a real helium receipt. It read the saved `-var DT 1e-7`
argument, selected step 40000 (4 ms), and produced all five figures. No
timestep override was supplied.

## Installed workflow

Follow the [README](../README.md) for a short solver installation check and
the 12 ms helium example followed by molecule tracing. A successful pipeline
has zero return codes and complete manifests for both stages. For the
default example, the converted field must identify timestep 120000 and the
molecule stage must retain exactly `N` birth and crossing records.

Mid-run analysis instead records the selected saved frame and the observed
helium status. A successful molecule calculation on that frame does not imply
that the parent helium calculation finished or reached equilibrium.

## Scope

Software checks establish that these inputs run and that the output pipeline
behaves as tested. They do not establish helium equilibrium, numerical
convergence, or agreement with an experiment. The 10 ms settling interval is
the tutorial starting prescription; molecular results use one frozen 2 ms
helium field. Quantitative studies need checks appropriate to their observables.
