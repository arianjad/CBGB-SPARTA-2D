# Sanity checks: cyl-3d ports scan (2026-10-05)

Hand calculations and small independent models checked against the BaF⁺ tracer
runs of the `ports` campaign. The test cases are 10, 15, 20, 30 and 40 SCCM, with
table-central cross sections and 2000 ions per point.

- Shared inputs: only the cross-section tables (`tracer/data/`) and the frozen
  DSMC neon field. Nothing is taken from the tracer code itself.
- Coordinates: x = 0 is the chamber side of the inlet wall. Ions start at
  x = −1.27 mm, the outer face of the inlet aperture.
- Code and data: scripts are in `tools/sanity/` (see its README for the awk
  commands that cut the density and flow fields from `frozen.grid`). Raw output is in
  `results/sanity/range_check.txt`, `range_check.png` and `more_checks.txt`.

## Summary

| Check | Result |
|---|---|
| Collisions to reach 1 eV | agree within ~1 collision |
| Stopping depth and entrance dE/dx | tracer bug found: ions go 2–4 mm (12–22%) too deep and the entrance dE/dx is ~25% low |
| Range scaling (depth × density) | constant within 5–9% across flows |
| Time to slow to 1 eV | agrees within ~1% |
| Thermal diffusion | theory mobility is 0.92 × the polarization limit; measured D within 0–23% at 20–40 SCCM |
| Extraction (drift–diffusion model) | within 1–4 points of the tracer |
| Neon flow | mass balance within 5%; chamber density lies between the free-molecular and continuum limits |

## Collision model used in the hand calculations

The tracer scatters isotropically in the centre-of-mass frame. BaF⁺ is 156.3 u and
Ne is 20.18 u, so γ = 4m₁m₂/(m₁+m₂)² = 0.405. One collision leaves the ion with the
energy fraction f = 1 − γu, with u uniform on [0, 1]. That gives:

- Mean energy loss ⟨ΔE⟩/E = 0.2025, so the stopping power is
  dE/ds = n σ(E) · 0.2025 · E.
- Typical (log-mean) ion ⟨ln f⟩ = −0.2373, so it takes
  N(E₀ → 1 eV) = ln(E₀/1 eV)/0.2373 collisions. The collision-to-collision spread
  gives about ±3.5 collisions.
- Direction: the mean forward velocity drops by m₁/M = 0.886 per collision and
  the speed by 0.888, so the direction cosine falls by 0.997 per collision. The
  projected depth is about 0.91 × the path length.
- Free path: 1/(nσ). The density n(x) is the on-axis DSMC density, from cells
  with r < 0.6 mm.

## Stopping depth, collisions and dE/dx: tracer bug found

| Case | Collisions to 1 eV, calc / sim | x at 1 eV (mm), calc / sim | dE/dx peak x (mm), calc / sim | Peak dE/dx (eV/mm), calc / sim |
|---|---|---|---|---|
| 10 SCCM, 0.5 keV | 26.2 / 26.9 | 20.3 / 23.2 | 1.0 / 5.5 | 51 / 36 |
| 15 SCCM, 2 keV | 32.0 / 32.8 | 23.0 / 26.4 | 1.4 / 4.2 | 158 / 123 |
| 20 SCCM, 1 keV | 29.1 / 29.9 | 13.1 / 15.2 | 0.6 / 1.7 | 136 / 104 |
| 30 SCCM, 5 keV | 35.9 / 36.6 | 17.0 / 20.7 | 1.4 / 4.0 | 456 / 354 |
| 40 SCCM, 3 keV | 33.7 / 34.7 | 10.0 / 12.0 | 0.7 / 2.5 | 457 / 340 |

What agrees: the collision counts, so the energy loss per collision is right. The
dE/dx peak is at the entrance in both (see "Why there is no Bragg peak" below). The
on-axis density is within 3% of the chamber average, so density is not the cause of
the mismatch.

Cause: `freePath` draws each flight with the gas density at the flight's start
point, then flies the whole distance (`tracer/ParticleTracing3D.jl`, the
`freePath(v, vrel, props[6], props[7])` call in `propagate`). Every ion starts in the
inlet bore, where the on-axis density is only 35–50% of the chamber value. At
10 SCCM it is 4.6e21 m⁻³ at x = −1.27 mm against 1.2e22 in the chamber. So the
first flight is about twice too long.

Confirmed with a standalone Monte Carlo of the same collision physics
(`tools/sanity/flight_mc.py`):

| Case | MC, tracer's rule: x at 1 eV / peak x | MC, exact sampling: x at 1 eV / peak x |
|---|---|---|
| 10 SCCM, 0.5 keV | 22.9 / 4.5 mm | 20.8 / 1.2 mm |
| 15 SCCM, 2 keV | 26.7 / 5.5 mm | 23.9 / 0.7 mm |
| 20 SCCM, 1 keV | 15.1 / 2.0 mm | 13.6 / 0.7 mm |
| 30 SCCM, 5 keV | 20.3 / 3.0 mm | 18.0 / 2.5 mm |
| 40 SCCM, 3 keV | 11.6 / 2.7 mm | 10.4 / 1.0 mm |

The tracer's rule reproduces the tracer; exact sampling reproduces the hand
calculation.

Effect: ions stop a few mm closer to the outlet than they should. Extraction is
probably somewhat overestimated. The bias also applies wherever flights cross
density gradients, such as the outlet channel. The old fluor-cell tracer uses the
same logic, so it is affected too.

Decision: keep the current ports campaign as is for consistency, and fix
future runs.

Fix: the new tracer option `--max-flight` caps each flight and redraws the free
path at the new point (a null step). This is exact because collisions are memoryless.
`tools/wsl/sbatch_cyl3d_tracer.sh` sets `MAXFLIGHT` to 0.2 mm (half the lattice)
for every run except RUNBASE `cyl-ports`, which keeps 0. A validation job
(`mf-test`, 20 SCCM, 1 keV, 400 ions, RUNBASE `sanity-mf`) is queued; it should give
x at 1 eV ≈ 13.1–13.6 mm and the dE/dx peak at ≈ 0.7 mm. To-do: make 0.2 mm the
tracer's own default after the campaign (see `docs/TODO.md`).

### Why there is no Bragg peak

A Bragg peak needs the ion to start above its stopping-power maximum, so that it loses
energy faster as it slows. Here the beam is far below both maxima:

- Electronic stopping for a Ba-like ion peaks around 1–5 MeV per nucleon, i.e.
  hundreds of MeV for BaF⁺. This is a rough scaling estimate.
- Nuclear stopping peaks at about 170 keV (reduced energy ε ≈ 0.3).
- The beam is 0.25–5 keV.

Below the maximum, both kinds of stopping fall as the ion slows:

- Nuclear: S = 0.2025 · E · σ(E), and the table's σ only grows about as
  E⁻⁰·²…E⁻⁰·⁵, so S ∝ E⁰·⁵…E⁰·⁸.
- Electronic: proportional to v, and only about 2% of the total here.

So dE/dx is largest at the entrance and falls with depth. This is not like a photon
beam: photons fall off exponentially because they are removed one at a time, each
keeping its full energy. Every ion here slows gradually, about 20% per collision, and
has a finite range. The tail is soft only because of straggling (the depth at 1 eV
spreads by ±4–7 mm). References: Lindhard, Scharff & Schiøtt (1963);
Ziegler, Biersack & Littmark (1985).

## 1. Range scaling

Depth from the beam start to 1 eV, times the chamber density, in 10²⁰ m⁻². At a fixed
energy it should be the same at every flow if only the density matters. Cells marked
"–" have fewer than 90% of ions reaching 1 eV.

| SCCM | 500 eV | 1000 eV | 2000 eV | 3000 eV | 5000 eV |
|---|---|---|---|---|---|
| 10 | 2.96 | 3.64 | – | – | – |
| 15 | 3.00 | 3.85 | 4.91 | – | – |
| 20 | 3.02 | 3.86 | 4.99 | 5.91 | 7.08 |
| 30 | 3.11 | 3.88 | 5.05 | 6.07 | 7.68 |
| 40 | 3.16 | 3.98 | 5.20 | 6.21 | 7.76 |

- Spread: 5–9% at each energy, rising slightly with flow, consistent with the
  free-path bias above.
- Monte Carlo with the tracer's rule: 2.97 (10 SCCM, 0.5 keV), 4.87
  (15 SCCM, 2 keV) and 3.83 (20 SCCM, 1 keV), matching the table.

## 2. Time to slow down to 1 eV

The Monte Carlo time is the sum of each free path divided by the ion's speed.

| Case | MC, tracer's rule (µs) | MC, exact (µs) | Tracer mean (µs) | Tracer median (µs) |
|---|---|---|---|---|
| 10 SCCM, 0.5 keV | 4.55 | 4.38 | 4.54 | 4.43 |
| 15 SCCM, 2 keV | 3.36 | 3.31 | 3.37 | 3.32 |
| 20 SCCM, 1 keV | 2.51 | 2.40 | 2.48 | 2.41 |
| 30 SCCM, 5 keV | 1.82 | 1.78 | 1.81 | 1.77 |
| 40 SCCM, 3 keV | 1.32 | 1.26 | 1.33 | 1.30 |

## 3. Thermal diffusion of BaF⁺ in 23 K neon

Theory, from the bare momentum-transfer table `bafplus_ne_q1.tsv` (Chapman–Enskog,
first approximation):

- Ω⁽¹'¹⁾(23 K) = 1.90e-18 m², which matches the table header.
- n·D = 2.56e19 m⁻¹ s⁻¹, so D = 1.09e-3 m²/s at 20 SCCM.
- The reduced mobility is K₀ = 4.80 cm²/Vs. The pure-polarization (Langevin) limit is
  K₀ = 13.876/√(α μ) = 5.22 cm²/Vs, so the ratio is 0.92. That is the expected size
  for 23 K.

Measured from the GIF trajectories (all energies of a flow, 40 ions each): the
sideways (y–z) mean-square displacement over 0.2–3 ms after t = 0.3 ms, using
MSD = 4 D Δt.

| SCCM | n_ch (m⁻³) | D theory (m²/s) | D measured (m²/s) | Ratio | Median ion x at 0.3 ms |
|---|---|---|---|---|---|
| 10 | 1.21e22 | 2.11e-3 | 5.8e-4 | 0.27 | 45.5 mm (at the end wall) |
| 15 | 1.77e22 | 1.44e-3 | 9.7e-4 | 0.67 | 30.0 mm |
| 20 | 2.34e22 | 1.09e-3 | 1.10e-3 | 1.01 | 23.8 mm |
| 30 | 3.50e22 | 7.3e-4 | 8.0e-4 | 1.10 | — |
| 40 | 4.69e22 | 5.4e-4 | 6.7e-4 | 1.23 | 10.6 mm |

- 20–40 SCCM: ions are mid-chamber and the measurement is clean. Agreement is
  within 0–23%; the remainder is flow-induced spreading.
- 10–15 SCCM: the ions sit against the end wall or in the flow converging on the
  outlet, which limits their sideways spread. So the low values there are not a clean
  diffusion test.

## 4. Extraction: independent drift–diffusion model

Model:

- Start: ions are slowed to 1 eV by the Monte Carlo with the tracer's flight rule,
  so they start from the same places as in the tracer.
- Motion: they drift with the DSMC gas velocity (averaged around the axis in
  0.4 mm bins) plus ∇D, and diffuse with the Chapman–Enskog D(n, T).
- Fate: walls absorb. An ion is "back" if it leaves through the inlet bore and
  "extracted" if it leaves through the outer face of the outlet channel.

| Case | Extracted, model / tracer | Back out the inlet, model / tracer | Wall, model / tracer | Median transit of extracted (ms), model / tracer |
|---|---|---|---|---|
| 10 SCCM, 0.5 keV | 3.3 / 5.3% | 1.1 / 0.6% | 95.6 / 94.0% | 11 / 18 |
| 15 SCCM, 2 keV | 11.3 / 14.8% | 0.8 / 0.1% | 87.9 / 85.0% | 17 / 14 |
| 20 SCCM, 1 keV | 10.9 / 13.5% | 10.5 / 3.8% | 78.6 / 82.8% | 39 / 36 |
| 30 SCCM, 5 keV | 32.4 / 33.1% | 5.4 / 1.7% | 62.2 / 65.2% | 35 / 29 |
| 40 SCCM, 3 keV | 24.0 / 27.8% | 40.0 / 16.0% | 36.0 / 56.3% | 47 / 43 |

- Extraction agrees within 1–4 points; the model runs slightly low. Transit times
  agree within about 20%.
- Inlet vs wall split at high flow: this is where the model and the tracer
  disagree. At 40 SCCM ions stop about 10 mm from the inlet, where the gas streams
  back out through the inlet at 20–120 m/s on the axis. The model's axis-averaged
  field smears out the one-sided gas inlet (on the −z side), so the tracer's 3D
  field is more trustworthy there.

Exit energy. The median exit energy of extracted ions is about 0.03 eV. That is
about 10 times 3/2·kT at 23 K, but it is the gas flow carrying them, not heat:

- Extracted ions are recorded at the edge of the simulation box (x = 78.4 mm).
- There they move forward at 160–185 m/s, against a gas speed of about 205 m/s. On the
  axis the gas speeds up from about 120 m/s at the outlet face to ~205 m/s; the limit
  for neon from 23 K is √(5kT/m) = 218 m/s.
- Their sideways speed is about 60 m/s, close to the BaF⁺ thermal speed at 23 K
  (56 m/s).
- So the ions are cold and trail the gas by about 15%. The 0.1 eV "thermalized" cut is
  unaffected.

## 5. Neon flow

Definitions:

- Q = inflow. Ratio = DSMC outflow/inflow at steady state.
- n_FM = 4Q/(c̄ A_eff): free-molecular, with the Clausing-weighted A_eff = 6.6 mm².
- n_cont = Q/(W c̄/4 (A_in + A_out)): ideal choked orifices with W = 1.82 × the
  free-molecular flux per area and no wall friction. This is the lowest density the
  flow could give.
- λ = Ne–Ne mean free path; Kn = λ / outlet diameter (3.175 mm).
- Re = outlet Reynolds number at the DSMC exit-plane density and speed.

| SCCM | Out/in | n DSMC (m⁻³) | DSMC / n_FM | DSMC / n_cont | λ (µm) | Kn | Re |
|---|---|---|---|---|---|---|---|
| 5 | 0.999 | 6.38e21 | 0.73 | 2.62 | 167 | 0.053 | – |
| 10 | 1.032 | 1.21e22 | 0.69 | 2.49 | 88 | 0.028 | 13 |
| 15 | 1.019 | 1.77e22 | 0.68 | 2.42 | 60 | 0.019 | 20 |
| 20 | 1.000 | 2.34e22 | 0.67 | 2.40 | 46 | 0.014 | 27 |
| 25 | 1.022 | 2.92e22 | 0.67 | 2.39 | 37 | 0.012 | – |
| 30 | 0.960 | 3.50e22 | 0.67 | 2.39 | 31 | 0.010 | 41 |
| 35 | 1.038 | 4.08e22 | 0.67 | 2.39 | 26 | 0.008 | – |
| 40 | 0.952 | 4.69e22 | 0.67 | 2.40 | 23 | 0.007 | 55 |
| 45 | 1.046 | 5.27e22 | 0.67 | 2.40 | 20 | 0.006 | – |
| 50 | 0.970 | 5.88e22 | 0.67 | 2.41 | 18 | 0.006 | – |

(Re uses an approximate viscosity of 2.7e-6 Pa s and is only given for the five flows
whose 2D field was extracted.)

- Mass balance: conserved within 5% at every flow.
- Density limits: the apertures are nearly continuum (Kn 0.006–0.05) but strongly
  viscous (Re 13–55), so wall friction is important. The DSMC density lies between
  the free-molecular and the frictionless-continuum limits, as it should.
- Not pinned down: the exact value needs transition-regime short-tube data, which
  this check does not include.
- Caveat: the apertures are only 3–4 lattice cells across at 0.4 mm, so the
  aperture conductance is the most lattice-sensitive quantity. It belongs in the
  lattice-convergence test in `docs/TODO.md`.
