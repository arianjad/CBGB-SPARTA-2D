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

## Portable checks

From the repository root after installation:

```bash
source .venv/bin/activate
bash tools/test.sh
```

These checks use synthetic data and temporary directories. They cover
clipped wall/domain contacts, exit statistics, first-crossing accounting,
trajectory coloring, geometry regeneration, and launcher refusal/failure
receipts. They do not download packages or launch SPARTA.

## Installed workflow

Follow the [README](../README.md) for a short solver installation check and
the 12 ms helium example followed by molecule tracing. A successful pipeline
has zero return codes and complete manifests for both stages. For the
default example, the converted field must identify timestep 120000 and the
molecule stage must retain exactly `N` birth and crossing records.

## Scope

Software checks establish that these inputs run and that the output pipeline
behaves as tested. They do not establish helium equilibrium, numerical
convergence, or agreement with an experiment. The 10 ms settling interval is
the tutorial starting prescription; molecular results use one frozen 2 ms
helium field. Quantitative studies need checks appropriate to their observables.
