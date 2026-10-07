# Handoff: tracing BaF+ through the completed 3D neon field

Follows `docs/hpc-fluor-cell-3d-handoff.md`. Input: a completed, settled
`cases/fluor-cell-3d` neon DSMC run (a `field.grid` with several complete,
non-transient frames). If that isn't done, finish it first; this doc doesn't
repeat the setup.

Question: does the hot (500 eV) BaF+ beam make it through the cell, how much
is lost, and how does that depend on buffer gas flow. Not a detailed
ion-collision physics study; see section 5 ("What this deliberately does NOT
model") before reading much into any one number.

## 0. What's new here

`tracer/ParticleTracing3D.jl` is a new file. `tracer/ParticleTracing.jl` stays
2D-axisymmetric only (used by `cases/b5-lean` and `cases/takahashi-single`)
and is untouched. The collision physics is reused unchanged: `collide!`, the
exact flux-weighted collision-partner sampler, and `freePath` are a general
elastic two-body treatment, valid for any mass ratio and relative velocity
(not a near-thermal approximation), so a 500 eV ion slowing in 23K neon needs
no new physics there. New: 3D geometry (segment vs triangle-mesh collision
detection against `cell_fluor3d.surf`, with a spatial grid so it isn't an
O(triangle count) scan per free-flight step) and a true 3D field lookup (no
assumed axisymmetry; the off-axis gas inlet breaks it).

Testing status: smoke-tested locally (single- and multi-threaded, a slow
near-thermal case and a fast 500 eV-equivalent case) on a synthetic cube
geometry and synthetic field, not on the real `cell_fluor3d.surf`/`field.grid`
(not available on the machine where it was written). Two bugs were found and
fixed in that testing (commented in the file where fixed): a zero-width
lookup-table range that crashed on near-uniform field data, and a
`Threads.threadid()` assumption that breaks with some Julia versions'
thread-pool models. The first real run is the first end-to-end test; go
through the section 4 checklist before trusting a result.

## 1. Freeze one field frame

`field.grid` appends one snapshot per `dump grid` interval over the whole
run, same as the 2D pipeline, which uses `tools/field2tracer.py` to pick one
frame for the tracer. The 3D tracer reads SPARTA's native column layout
directly, so the coordinate remapping and He-specific mass-density conversion
aren't needed and a simpler tool does it:

```bash
python tools/field3d_freeze.py "$DSMC_RUN_ROOT/<your-run-name>" \
  --out cases/fluor-cell-3d/frozen_field.grid
```

Default is the latest complete frame (same selection flags as
`field2tracer.py`: `--timestep`, `--until-step`, `--until-ms` + `--dt`).
Make sure the frame isn't still transient (README: "Allow the first 10ms for
settling").

## 2. Smoke test against the real geometry

Run a handful of particles first and look at the output:

```bash
cd cases/fluor-cell-3d
julia --project=../../tracer --threads=4 ../../tracer/ParticleTracing3D.jl \
  cell_fluor3d.surf frozen_field.grid \
  -n 20 -x 0.03048 -y 0.01905 -z 0.01905 \
  --vx 24843.71 --vy 0 --vz 0 \
  -m 20.1797 -M 156.325 \
  --saveall 1 --seed 1
```

(`--sigma-model table` is the default: the computed BaF+-Ne cross-section
curve in `tracer/data/bafplus_ne_sigma_tracer.tsv`, see section 3. The older
models are `--sigma-model edep` and `--sigma-model const --sigma <value>`.)

Check:
- It runs without errors. A range/bounds crash is a real-data edge case the
  local tests didn't hit; the comments at the two fixed bugs in
  `ParticleTracing3D.jl` show what that kind of bug looks like. Don't just
  silence the error.
- `--spawnclip` (on by default) didn't reject the spawn point as inside
  geometry. (0.03048, 0.01905, 0.01905) is the molecule-inlet port location
  from `gen_fluor3d.py`. If the geometry changed and the port moved, every
  particle fails to spawn and it errors out after 1000 tries.
- Exit positions and velocities make physical sense (section 4), not just
  "it produced numbers".

## 3. Physics parameters

| Parameter | Value | Why |
|---|---|---|
| `-M` (BaF+ mass) | `156.325` amu | Matches the existing `cases/b5-lean/baf-hot-source.conf` (156.325 u) |
| `-m` (buffer gas mass) | `20.1797` amu | Neon, matches `ne.species` |
| `--vx` | `24843.71` m/s | v = sqrt(2 * 500eV / m_BaF+); recompute if the injection energy changes |
| `-x -y -z` | `0.03048 0.01905 0.01905` | Molecule-inlet port location (2.54mm channel, axial) |
| `--sigma-model` | `table` (default) | Computed BaF+-Ne cross section vs collision energy, see below |
| `--sigma-table-variant` | `central` (default) | `low`/`high` bracket the short-range-wall uncertainty (only matters above ~1 eV) |
| `-T` | `0.0` (default) | Monoenergetic, mono-directional beam. A simplification, not a measurement: if the real source has a known energy/angular spread, set a nonzero `-T` (thermal-style spread around the mean velocity) or check whether the spread matters for the conclusion. |

### Cross section, current (supersedes the older notes below)

The default is `--sigma-model table`, a computed BaF+-Ne cross section over
the whole 500 eV-to-thermal range (`tracer/data/bafplus_ne_sigma_tracer.tsv`;
derivation and scripts in `tools/cross_section/`):

- Physics: classical momentum-transfer cross section Q1(E) for a Ba+-Ne proxy
  potential: the CCSD(T)/CBS Ba+-Ne potential of Buchachenko & Viehland,
  J. Chem. Phys. 148, 154304 (2018) (well + exact ion-induced-dipole C4 = alpha_Ne/2)
  joined to the ZBL universal repulsion at short range. It reproduces the
  1.105 x Langevin polarization limit and agrees with ZBL/SRIM universal
  nuclear stopping to 2-15% at 100-500 eV lab.
- Tracer mapping: the file is not bare Q1. The tracer uses a flux-weighted
  partner + isotropic scattering at rate n*sigma*v_coll, so the table is the
  drag-matched sigma_t(E_cm) that makes the tracer's mean drag equal the true
  drag at every ion speed. It equals Q1 for fast ions and gives the correct
  thermal diffusion coefficient at rest (2.0e-18 m^2 at 23 K). Bare Q1 is in
  `tracer/data/bafplus_ne_q1.tsv` for reference only; don't pass it to `--sigma-table`.
- Key values (sigma_t): 2.0e-18 m^2 thermal (23 K), 1.0e-18 at E_cm 100 K,
  2.2e-19 at 1 eV lab, 1.4e-19 at ~8 eV lab, 4.9e-20 at 500 eV lab.
- Synthetic check (uniform 23 K, 2e22 m^-3 Ne, 500 eV injection, 20 ions): mean
  stopping depth 14 mm (table) vs 8.8 mm (old edep), 11 mm (const 1.5e-19),
  0.8 mm (const 2e-18); ~40 collisions to thermalize.
- Uncertainty: thermal/low-energy end +-10-15% (set by the long-range
  physics); above ~1 eV +-30-40% (short-range wall), bracketed by
  `--sigma-table-variant low` and `high`. No BaF+-Ne cross section has been
  measured; Ba+ stands in for BaF+.
- The default table is only valid for BaF+ in Ne; the tracer refuses to start
  with it for other masses.

Run the `central` table as the main condition, plus `--sigma-table-variant low`
and `high` as the bracket. The older notes below are kept for history.

### Cross section, earlier revision (relevant to transmission numbers from older versions of this doc)

The first version used a single constant `--sigma`, scanned between two
unrelated guesses (5.0e-19 and 1.5e-19 m^2); transmission changed a lot
between them (72% vs. 89%), an open uncertainty at the time. Two changes
since:

1. Energy-dependent cross section (`--sigma-energy-dependent`, on by default)
   instead of a constant: `sigma-low` at/below the local buffer-gas thermal
   speed, falling as 1/v (Langevin capture: for a polarization-dominated
   ion-neutral potential the capture cross section goes as 1/v, so the rate
   constant sigma*v is velocity-independent) down to a `sigma-high` floor at
   high relative speed. The same two numbers become the endpoints of one
   curve rather than an arbitrary bracket; see `sigma_of_speed()`'s docstring
   for the derivation. The constant mode is still there
   (`--sigma-energy-dependent 0 --sigma <value>`) for comparison.
2. `sigma-low` raised from 5.0e-19 to 2.0e-18 m^2 (4x), from measured data:
   two buffer-gas-cooling papers measured the BaF (neutral) + He collision
   cross section directly, 1.4e-18 m^2 (*Phys. Rev. A* 95, 032701 (2017)) and
   2.7e-18 m^2 (arXiv 1906.08798), in good agreement. The old 5.0e-19 was an
   indirect extrapolation (Ba+ ion mobility in He/Ar, not BaF, converted
   through a mobility formula and then polarizability-scaled to Ne). The new
   value is scaled from a direct measurement of the actual molecule in He, so
   one scaling step (polarizability, He->Ne) instead of two. It's still the
   neutral molecule, not the BaF+ ion; ion-induced-dipole attraction should
   make the ionic value equal or larger, so 2.0e-18 is closer to a lower
   bound than a central estimate.

`sigma-high` (1.5e-19 m^2) is unchanged and still just an estimate from
summed ionic/atomic radii. No eV-keV data for this or a close analog
ion-neutral system turned up in either literature search. One loose hint from
general ion-atom scattering (a shallower E^(-1/3) power law reported for the
elastic-scattering amplitude in some systems, not confirmed here) suggests the
high-energy cross section may fall off more slowly than 1/v, i.e. stay higher
over more of the trajectory. Not a number to plug in, just a reason not to
over-trust the high-energy end either.

Since both ends are uncertain, run and report at least three conditions: the
default energy-dependent model and the constant-sigma pair
(`--sigma-energy-dependent 0 --sigma 2.0e-18` and
`--sigma-energy-dependent 0 --sigma 1.5e-19`) as bounds. If the flow-rate
conclusion is the same for all three, the cross-section uncertainty doesn't
block the decision; if not, that matters.

### Output: success and loss

Columns: `idx x y z xnext ynext znext vx vy vz collides time`;
`xnext,ynext,znext` is where the path ended.
- Transmitted: `xnext` near 0.08255 m (molecule-outlet port, `MOLECULE_OUTLET`
  in `gen_fluor3d.py`) and `collides` usually large (it thermalized and
  diffused out rather than coasting straight through; worth checking, see
  below).
- Lost back out the inlet: `xnext` near 0.03048 m again.
- Stuck on a wall: `xnext` on the housing wall (not at either port's x/y/z).
  Same convention as 2D, molecules stick when they hit a wall:
  `collide!`/`propagate` returns as soon as `getCollision` reports a wall
  hit, so the final `xnext,ynext,znext` is where it stuck.

Transmission fraction = (# ending near the outlet) / `n`. `--saveall 1`
prints every particle regardless of outcome (needed for that ratio). The
default (0) prints only particles that left the simulation box, a box-edge
criterion that undercounts "reached the outlet" when the box extends past it
(it does, with the scoping deck's placeholder `create_box`).

## 4. Validation checklist (before trusting a real result)

1. Energy: for a sample of transmitted/lost particles, the exit speed
   (`sqrt(vx^2+vy^2+vz^2)`) should be far below the 24843.71 m/s injection
   speed, comparable to the neon field's own velocity scale (tens to
   hundreds of m/s). Particles exiting near injection speed with very few
   `collides` mean something is wrong (density/cross-section scale, or the
   particle punched through without seeing the gas); don't use a
   transmission fraction from that run.
2. Collision count: compare `collides` per particle with a rough estimate:
   mean free path lambda = 1/(nrho*sigma) with the buffer gas density near
   the inlet from the neon run, and the cell scale (~50mm). An order of
   magnitude off from a simple diffusion estimate
   (steps ~ (distance/lambda)^2) needs looking into.
3. Geometry read: confirm `read_surf3d` parsed the expected point/triangle
   counts, e.g. a quick `println(length(points), " ", length(triangles))`
   after the `read_surf3d` call, compared with `check_surf3d.py`'s counts
   for the same file.
4. Reproducibility: same `--seed`, same thread count -> same result (as in
   the 2D file, threaded runs only reproduce at a fixed thread count, not
   across `--threads` values). If a run doesn't reproduce from its recorded
   seed, don't trust it.
5. Run all three cross-section conditions from section 3 (energy-dependent
   default and both constant-sigma bounds) and report all three.

## 5. What this deliberately does NOT model

From the discussion that scoped this work, restated because it's easy to
forget once the pipeline produces numbers:

- No dissociation / internal energy exchange. BaF+ is a single rigid point
  mass. One hard early collision can plausibly transfer more than BaF's
  ~5-6 eV bond dissociation energy; the model can't represent breakup and
  assumes the molecule survives the whole trajectory. It can't answer
  survival-fraction questions, only transport.
- No charge exchange or other reactive channels.
- Energy-independent cross section (one `--sigma` over the whole
  500eV-to-thermal slowing-down range) rather than a modeled
  energy-dependent one. The section 3 sensitivity scan mitigates this but
  doesn't fix it.
- No ion guide / applied field. Free-flight legs are straight lines (the 2D
  file's harmonic-trap-aware propagation doesn't apply since this cell has no
  such field; confirm that's true for the actual setup before trusting a
  trajectory shape).
- Frozen buffer-gas field (same one-way coupling as the 2D pipeline). A
  500 eV ion's energy deposition into the gas is assumed negligible against
  the gas heat capacity at the actual ion flux; not checked against a real
  ion rate.
- `StatsArray`'s r/z binning is a 2D summary projection (axial x vs.
  distance off the x-axis), kept for continuity with the 2D output. It won't
  show y/z asymmetry from the off-axis gas inlet; fine for a quick profile,
  not for anything needing 3D resolution.

## 6. What to report back

- Whether the smoke test (section 2) and a statistics-quality run (larger
  `-n`, at least a few hundred particles per condition) completed cleanly.
- Transmission / lost-to-inlet / stuck-on-wall fractions for all three
  cross-section conditions from section 3, at the flow rate(s) of the
  completed neon run(s).
- Section 4 checklist results: whether the numbers hold together
  physically, not just "it ran".
- For a flow scan: which neon runs (or new ones at different
  `MDOT`/`FILLN`) were traced against, and the transmission vs flow trend.
