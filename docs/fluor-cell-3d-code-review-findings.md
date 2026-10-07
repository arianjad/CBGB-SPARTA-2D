# fluor-cell-3d code review findings

Correctness check of `cases/fluor-cell-3d/` (SPARTA/neon side) and
`tracer/ParticleTracing.jl` (BaF+ ion side) before spending HPC time on the
buffer-gas run. Repo state: branch `experiments`, commit `511a0bc` (after the
duplicate-edges fix, the part-fusing fix, and the stray-optical-window fix).

## A. Neon/SPARTA side

Checked, OK:

- `MDOT`/`FILLN` aperture-law derivation in
  `docs/hpc-fluor-cell-3d-handoff.md` step 4. Recomputed from scratch
  (v_mean(23K Ne) = 155.36 m/s, MDOT(4 sccm) = 6.0022e-8 kg/s,
  FILLN ≈ 5.44e21); matches the doc's values exactly. Same method as
  `cases/b5-lean/gen_b5.py`'s He derivation, redone for Ne.
- Ne atomic mass in `ne.species` (3.350918e-26 kg) = 20.1797 amu ×
  1.66053906660e-27 kg/amu.
- `mixture negas Ne nrho 1.0e22 temp 23.0` in `in.ne_fluor3d_scoping` looks
  like a placeholder but follows the same convention as
  `cases/b5-lean/in.he_b5_mflow` and `cases/takahashi-3d/in.he_takahashi_3d_scoping`.
  The real flux comes from `emit/surf ... mflow`, not this nrho. Intentional,
  not a copy-paste bug.
- Geometry bounding box sits inside the `create_box`/region `gas` bounds in
  both the smoke-test and scoping decks, with reasonable margin.
- `gen_fluor3d.py`'s mm to m conversion (`MM_TO_M = 1e-3`) is applied once, at
  final point conversion. No unit mixing.
- `rotdof=0` in `ne.species` is right for a monatomic gas (dof=3,
  translational only), consistent with Bird 1994 Table A1's dof=3 entry for Ne.

Bug (small, quantified):

- `ne.vss` VHS reference diameter is off by ~1%. Recomputing the VHS
  diameter from Bird 1994 Table A1's own inputs for Ne (mass = 3.350918e-26
  kg, reference viscosity μ_ref = 2.975e-5 N·s/m² at T_ref = 273.15 K,
  omega = 0.66) with the standard VHS relation

  ```
  mu_ref = 15 * sqrt(pi * m * kB * T_ref) / (2*pi * d_ref^2 * (5-2*omega)*(7-2*omega))
  ```

  gives d_ref ≈ 2.766e-10 m, not the `2.79e-10` in
  `cases/fluor-cell-3d/ne.vss`. That's ~0.9% in diameter, ~1.7% in
  cross-section (sigma ∝ d²). Small, but it follows directly from the table's
  own numbers, not rounding. Action: set the `ne.vss` diameter to `2.766e-10`
  (T_ref = 273 vs 273.15 K only changes the 4th significant figure, so that
  isn't the source).

Already known, now quantified:

- The flat `create_grid 60 27 33` in `in.ne_fluor3d_scoping` gives cells
  ≈3.25–3.33mm on a side, vs. the smallest port feature (1.588mm gas-inlet
  cap), so at most 1–2 cells across that port. The handoff doc (step 5) and
  the deck's own TODO comment already flag this grid refinement gap. It is a
  real resolution problem, but per the doc it doesn't block a first run.

## B. BaF+/Julia tracer side: what isn't built yet

1. `tracer/ParticleTracing.jl` can't read this case's output.
   `build_field` parses a 2D-axisymmetric 8-column dump layout
   (`x,y,T,rho,rho_m,vx,vy,vz`) matching the He cases'
   `dump gd grid all ... id xc yc xlo ylo xhi yhi f_ag[*]`.
   `interpolate!` rebuilds a 3D velocity from one stored radial component via
   `interp[4]*cos(theta), interp[4]*sin(theta)`, i.e. assumes azimuthal
   symmetry. `in.ne_fluor3d_scoping`'s dump
   (`id xc yc zc xlo ylo zlo xhi yhi zhi f_ag[*]`, full 3D u,v,w) has a
   different column layout, and the flow is non-axisymmetric, which is the
   reason this case exists. `getCollision` only intersects 2D r-z line
   segments (2D `.surf` format), not the 3D triangle mesh. A 3D field
   reader/interpolator/collision detector doesn't exist yet (the handoff doc
   says so too). This is a hard blocker.

2. The elastic hard-sphere collision model isn't valid at 500 eV, regardless
   of dimensionality. `collide!` / `exact_partner` do energy-conserving
   elastic scattering with a fixed cross section (CLI scalar `--sigma`): no
   energy dependence, no polarization/Langevin term, no charge exchange.
   `in.ne_fluor3d_scoping`'s header already says this: a 500 eV BaF+ ion
   thermalizing in 23K neon is a stopping-power problem (many collisions,
   large fractional energy loss each), "flagged, not implemented." Spawning a
   monoenergetic 500 eV particle already works (`-T 0 --vz <speed>` gives a
   delta-function velocity); the collision physics after that is missing.

3. `cases/b5-lean/baf-hot-source.conf` is not a usable template.
   `SOURCE_T_K=1000.0` (thermal ablation source, ~0.1 eV mean energy, four
   orders of magnitude below 500 eV), different buffer gas (He, not Ne),
   different cell. Its `CROSS_SECTION_M2=2.7e-18` was tuned for BaF-He at
   ~1000K and must not be reused for BaF+-Ne at 500 eV. Flagged because it's
   the obvious file to copy from.

4. No charge/Coulomb code exists anywhere in the tracer.
   `ne.species`'s charge column is 0.0 (correct, SPARTA only sees neutral
   Ne). A polarization or charge-exchange treatment for the ion side would be
   new code.

Summary: the neon/SPARTA side is fine after the VHS diameter fix; everything
else checked out numerically, with one already-tracked grid-resolution gap
that doesn't block a first run. The BaF+ side needs (a) a new 3D field
reader/interpolator/collision detector and (b) a real stopping-power
collision model before any 500 eV ion trajectory means anything. Plotting the
neon field works now; BaF+ results are a separate follow-up.
