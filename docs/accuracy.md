# Controlling numerical accuracy

There is no universal 5% accuracy setting. Choose an observable and a tolerance:
for example, 1% relative uncertainty in extraction, or 1 m/s in mean axial speed.
The tutorial demonstrates a working calculation; it does not establish either
5% or 1% accuracy for every output. These controls are available for exploration;
the convergence study is needed when making a quantitative claim.

## What to change

| Control | Where to change it | What the comparison tests |
| --- | --- | --- |
| Helium grid | `create_grid` and refinement `region` definitions in `cases/b5-lean/in.he_b5_mflow` | Spatial resolution of collisions, gradients, inlet and aperture flow |
| Helium simulation particles | `FNUM` in the helium launch | Finite-particle noise and particle-count dependence of mean fields |
| Helium timestep | `DT`, with corresponding `STEPS` and averaging settings | Time discretization at fixed physical duration |
| Settling and sampling | `STEPS`, `fix ag ave/grid`, station averages, saved-field selection, independent helium `SEED` values | Residual transients, temporal variation, and uncertainty of averaged quantities |
| Molecule sample size | `N` and `SEED` in `tools/run_2d_standalone.sh` | Molecular Monte Carlo uncertainty conditional on the selected helium field |

`MDOT`, wall temperature, source parameters, and cross sections describe the
physical problem. Vary them to study sensitivity or input uncertainty, keeping
them fixed during a numerical-convergence comparison. More MPI ranks or Julia
threads primarily change cost and random-number execution, not resolution.

## Helium particle count: lower FNUM

At fixed grid, weighting, density, and flow, reducing `FNUM` increases the
number of simulated helium particles. For example, `2.5e17` to `6.25e16`
requests approximately four times as many. The B5 deck uses radial cell
weighting, so a particle's physical weight also depends on its cell; `FNUM`
alone is not its universal physical multiplicity. Inspect actual particle
counts, including sparsely populated refined cells. See SPARTA's
[global fnum and weighting documentation](https://sparta.github.io/doc/global.html).

In the public launcher, a particle-count comparison can start with:

```bash
FNUM=2.5e17 SEED=101 bash tools/run_helium.sh particles-base
FNUM=6.25e16 SEED=101 bash tools/run_helium.sh particles-4x
```

Repeat with independent helium seeds. Four times as many independent samples
would halve sampling error, but temporal correlations and weighting mean this
is a planning estimate, not a promised reduction. Check both scatter and mean
changes. The same seed across different meshes or particle counts does not
produce paired trajectories. Increasing molecule `N` cannot repair a noisy
helium field.

## Mesh refinement and adaptive grids

The supplied mesh is **static with local refinement**, not an automatically
adapting mesh. Its `create_grid 244 50 1 levels 5 ...` and nested regions set
the resolution; each `2 2 1` subdivision halves both in-plane cell dimensions.
Refine or widen regions around the inlet, aperture, and steep gradients, then
compare the same physical observables. A finer grid needs adequate particles
in its smaller cells and may need a smaller timestep. See
[create_grid](https://sparta.github.io/doc/create_grid.html).

Keep the physical geometry fixed. The geometry generator's `DX`, `DR`, and
refinement-level checks must match the deck; its cap coordinate `XC` is itself
derived from `DX`. Blindly changing `DX` and regenerating would also move a
physical surface. Preserve wall locations and inspect wall/station alignment
when redesigning the grid. Keep the nested transition bands and use
`tools/check_grid_2d.py` on a frozen snapshot to inspect the current 2:1
adjacency convention. That check inspects one frame and does not establish
mean-free-path resolution or accuracy.

For physically informed refinement, SPARTA's
[compute lambda/grid](https://sparta.github.io/doc/compute_lambda_grid.html)
provides mean free path, mean collision time, and directional ratios
`knx = lambda/dx`, `kny = lambda/dy`. Inspect populated regions relevant to
the observable; small cell Knudsen numbers identify cells large relative to
the mean free path. Assess gradients and geometry as well. No single threshold
guarantees a percentage error in extraction.

To record these diagnostics, insert the following after the existing
`fix ag`/`dump gd` definitions and before `run`:

```text
compute resolution lambda/grid f_ag[1] f_ag[4] lambda tau knx kny
dump resolution grid all 20000 resolution.grid id xc yc c_resolution[*]
```

The four diagnostic columns are `lambda` (m), `tau` (s), `knx`, and `kny`.
Keep them in this separate dump so the converter's `field.grid` columns remain
unchanged. Match the dump period to `ag` if changing its frequency. Interpret
populated averaging windows, not the empty initial frame. This snippet passed
a syntax/output check with the pinned solver; it does not set an accuracy target.

SPARTA also supports a one-time [adapt_grid](https://sparta.github.io/doc/adapt_grid.html)
operation and periodic [fix adapt](https://sparta.github.io/doc/fix_adapt.html).
They can refine/coarsen by particle count, surface proximity, or a computed
per-cell quantity. For a Knudsen criterion, smaller `lambda/cell-size` calls
for refinement (`thresh less more`); choose thresholds and `maxlevel` from
resolution diagnostics and cost. Particle-count refinement alone is not a
mean-free-path criterion.

Automatic adaptation requires a custom deck; there is no tested automatic
accuracy mode in this example. A practical route is to adapt during a pilot,
retain the resulting grid, and collect statistics on a fixed mesh. With
on-the-fly adaptation, averaging windows must not span mesh changes; SPARTA
restricts `fix ave/grid` accordingly. The included plotter can inspect a
single complete adaptive-mesh frame, but it refuses to average frames with
different cell IDs rather than inventing a remapping. See
[ave/grid restrictions](https://sparta.github.io/doc/fix_ave_grid.html).

## Timestep and physical averaging windows

Compare `DT` and `DT/2` at equal physical duration. For the tutorial,
`DT=5e-8 STEPS=240000` still means 12 ms. Cell transit times and collision
times both matter: SPARTA's
[compute dt/grid](https://sparta.github.io/doc/compute_dt_grid.html)
estimates candidate timesteps from both, with user-selected fractions.
Use those diagnostics and a timestep-refinement comparison, rather than
assuming 100 ns is sufficient for a changed flow or mesh.

Also preserve the averaging window when comparing timesteps. The supplied
`fix ag ave/grid all 10 2000 20000 ...` samples every 1 microsecond and saves
every 2 ms at `DT=1e-7`. Halving `DT` without editing it changes those times
to 0.5 microseconds and 1 ms. To preserve both times, change its three integers
to `20 2000 40000` and change `dump gd ... 20000` to `40000`. Scale station
averages/dumps similarly if comparing their time series. See
[ave/grid sampling rules](https://sparta.github.io/doc/fix_ave_grid.html).

The time selectors read the saved fixed `DT`; they do not reconstruct a
variable-timestep history. If adding `fix dt/reset` or timestep resets to a
custom deck, record actual elapsed times and adapt the analysis before using
millisecond selection or labels. Exact step selection remains available.

## Molecule count and helium-field uncertainty

For independent molecular trials on one fixed helium field, an extraction
probability `p` has binomial standard error `sqrt(p*(1-p)/N)`. Its relative
standard error is `sqrt((1-p)/(N*p))`. At `p=0.05`, the planning counts are:

| Target relative standard error | Approximate N |
| --- | ---: |
| 5% | 7,600 |
| 2% | 47,500 |
| 1% | 190,000 |

These are one-standard-error targets for that probability, not 95% confidence
bounds or guarantees of total accuracy. At `N=1000`, a 5% yield has about
14% relative sampling error (0.69 percentage points absolute). Rare accepted
subsets need larger samples; speed widths and tail probabilities need their
own uncertainty estimates.

Use the same source and acceptance definition on several helium snapshots
and independent helium runs. Separate molecule-draw noise from field-to-field
variation. Adjacent snapshots may be correlated; use time blocks long enough
to assess that correlation. Longer averaging should follow settling, not
hide drift. Tracing through a mean helium field is not generally equivalent
to averaging results traced through individual fields, because transport is
nonlinear in the field.

## A practical convergence sequence

1. Choose the observable, physical sampling window, and relative or absolute
   tolerance. Save the baseline configuration and uncertainty estimate.
2. Increase helium particles, using multiple seeds, until statistical scatter
   is small enough to distinguish changes at the requested tolerance.
3. Refine the mesh and reduce the timestep separately, then check the combined
   refined settings. Keep geometry, flow, duration, and averaging windows fixed;
   maintain adequate cell populations as the mesh changes.
4. Increase molecule `N` and compare independent helium fields. Require changes
   between successive refinements to be small relative to the target, with
   uncertainties small enough to make that comparison informative. More than
   one refinement level helps avoid an accidental agreement.
5. Report the remaining model/input uncertainties separately: gas collision
   parameters, wall interaction, molecular source and cross section,
   axisymmetry, frozen-field coupling, and the tracer's mean-speed collision
   rate and per-flight field sampling. Mesh refinement and larger `N` do not
   remove these approximations.

These are user-controlled studies. The launchers do not demand a particular
grid, `FNUM`, runtime, molecule count, or convergence verdict before allowing
an exploratory calculation.
