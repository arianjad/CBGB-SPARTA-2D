# fluor-cell-3d BaF+ tracer: bugs found and validated results

> **SUPERSEDED (2026-09-23): the neon field these results were traced
> against was wrong.** The surf file had the optical-window slot in the
> housing's top face left open, and ~70% of the injected neon escaped through
> it. The chamber sat at ~2e20 m^-3 instead of ~5.6e21. Everything below the
> "Sealed-window rerun" section was computed on that leaky field and is kept
> only for the record. See [Sealed-window rerun](#sealed-window-rerun) for
> current results.

## Sealed-window rerun

### How the leak was found

The exterior plume was lopsided upward. Measured over half-spaces about the
port axis (y = z = 19.05mm), in a |dy|,|dz| < 38mm window:

| Region | density +z/-z | x-flux +z/-z | control: density +y/-y | control: flux +y/-y | density-weighted mean w |
|---|---|---|---|---|---|
| 2.5-26.5mm past outlet | 4.25 | 5.02 | 0.90 | 0.86 | +31.9 m/s |
| 26.5-117.5mm past outlet | 2.00 | 2.46 | 0.92 | 0.91 | +9.4 m/s |
| behind the back port | 3.10 | 3.86 | 0.82 | 0.81 | +28.7 m/s |

**Ruled out: inlet placement.** Two controlled Takahashi-cell 3D runs gave a
symmetric plume (+z/-z within ~1-5%, the same as their +y/-y controls):

- `cases/takahashi-3d-offset`: inlet tube axis shifted 3mm off the cell axis.
- `cases/takahashi-3d-side`: inlet tube entering through the side wall from
  -z, like this cell's gas inlet.

**Found: an open slot in the top face.** `gen_fluor3d.py` excludes the
optical window (#1/#4) on the assumption that the top face is continuous. It
isn't. A 30 x 7mm slot (x 42-72mm, y 15.4-22.6mm, 229mm^2) runs from z=38.1
down to the chamber ceiling (z~34.3). Evidence of the leak:

- Gas rose inside the chamber, from ~0 m/s at the floor to +16 m/s just
  under the slot.
- It left the slot at ~+82 m/s.
- The upward flux through a plane just above the housing was 1.3e18 /s,
  against 1.79e18 /s injected (4 SCCM).

A matching bottom slot is sealed by the z=-6.35 plate, which is why only the
top leaked.

### Fix

`cases/fluor-cell-3d/seal_top_slot.py` edits `gen_fluor3d.py`'s output
directly (no cadquery needed) and writes `cell_fluor3d_sealed.surf`:

- It removes the slot's 24 side-wall triangles.
- It caps the two rims left behind (the top-face rim and the
  chamber-ceiling rim) with fans wound to match their neighbours.
- The result has 0 free, 0 non-manifold and 0 inconsistently wound edges.
  The solid volume grows by 950mm^3, which is the slot's volume.

Sealed, the chamber fills to the handoff's aperture-law density. There, the
23K Ne mean free path is ~0.4mm, so the flat 3.25mm grid would be ~9 mean
free paths per cell. `in.ne_fluor3d_sealed_run` therefore refines the
chamber and port channels 8x8x8, to ~0.41mm cells (1.49M cells total). Other
run settings:

- FNUM 2e11
- DT 5e-7 (still < 1/4 of the mean collision time)
- 50000-step averaging windows

The run took 600000 steps (0.3 s), as six ~9-minute mit_quicktest jobs
(`tools/wsl/sbatch_fluor3d_sealed.sh`).

### Sealed neon field (averaged over steps 450000-600000)

- **Chamber:** <n> = 5.63e21 m^-3, uniform to <1% with height (aperture-law
  estimate: 5.4e21); T = 23.0 K; all mean velocities < 0.3 m/s.
- **Settling:** total nrho drifts < 1% over the last 250000 steps.
- **Above the housing top:** 4e14 m^-3, down from 3.5e18.
- **Averaging caveat:** a resumed job re-dumps its starting step with an
  empty average. Drop duplicate timesteps before tail-averaging; leaving them
  in lowered averages by ~20%.

| Region | density +z/-z | x-flux +z/-z | control: density +y/-y | control: flux +y/-y | mean w |
|---|---|---|---|---|---|
| 2.5-26.5mm past outlet | 0.833 | 1.243 | 0.830 | 1.247 | -0.36 m/s |
| 26.5-117.5mm past outlet | 0.916 | 1.094 | 0.910 | 1.101 | -0.96 m/s |
| behind the back port | 0.820 | 1.271 | 0.820 | 1.277 | -0.22 m/s |

The plume is now symmetric: z matches the y control to 3 decimals. The ratios
aren't 1.00 because the exterior grid (3.33mm) doesn't split at the axis, so
the cells straddling the axis are all counted in the minus half. The spacing
is the same in y and z, so the bias is the same in both columns.

### BaF+ tracer on the sealed field

Setup: 1000 particles per condition, `--seed 42`, frozen step-600000 frame,
500 eV (24843.71 m/s) injected at the molecule-inlet port. The three
conditions from handoff section 3 (commit 5bca552):

| Cross section | Transmitted | Lost (back out inlet) | Wall hit | Mean collisions |
|---|---|---|---|---|
| energy-dependent (default: 2.0e-18 falling as 1/v to a 1.5e-19 floor) | **0.8%** (8) | 0.1% | 99.1% | 10640 |
| constant 2.0e-18 | 0.0% | 35.1% | 64.9% | 3924 |
| constant 1.5e-19 | 0.3% (3) | 0.0% | 99.7% | 50 |

In the leaky field, the 1.5e-19 and 5.0e-19 runs gave 72-89% transmitted.

At the real chamber density, a 500 eV BaF+ thermalizes within a few mm of
entry. It then diffuses to a wall (the tracer treats a wall hit as a loss)
long before the near-stagnant gas (< 0.3 m/s) can carry it to the 3.175mm
outlet. By condition:

- **Energy-dependent:** wall hits are spread along x = 48-77mm, with a pile
  on the end wall around the outlet channel at x = 77.2mm.
- **Constant 2.0e-18:** ions stop right at the entrance (wall hits at
  x ~ 31-33mm), and 35% diffuse back out the inlet channel.
- **Constant 1.5e-19:** ions coast most of the way and then strike the end
  wall around the outlet.

Thermalization check: ions with > 1000 collisions leave at a mean 74 m/s
(energy-dependent) or 77-87 m/s (constant 2e-18). The flux-weighted mean
speed at 23 K is 66 m/s, so there is no heating artifact. The higher overall
mean exit speed in the energy-dependent run (232 m/s) comes from ions that
hit a wall before thermalizing.

Scripts used (classifier, half-space asymmetry analysis) are in
`tools/fluor3d_asymmetry.py` and `tools/tracer_classify3d.py`.

## Flow-rate scan, sealed cell (2026-09-23 to 26)

### Setup

**Neon.** 18 flows from 1 to 58 SCCM, each an independent sealed-cell run of
0.3 s simulated.
- Same deck as 4 SCCM (`in.ne_fluor3d_sealed_run`, 8x8x8-refined chamber at
  0.41mm).
- MDOT, the aperture-law initial fill and FNUM all scale with SCCM, so each
  run has ~1M particles.
- DT is <= ~0.5 x the Ne-Ne collision time. Per-flow values are in
  `tools/wsl/fluor3d_flow_params.sh`.
- Every flow settled: 23 K, a symmetric plume, and a chamber density that
  rises almost linearly with flow (1.6e21 at 1 SCCM to 7.3e22 m^-3 at
  58 SCCM).
- 127 SCCM was also run but is excluded: ions don't leave via the outlet
  there.
- **From 10 SCCM up the chamber grid is coarser than twice the Ne mean free
  path**, so the port flow there is approximate.

**Tracer.** 500 eV BaF+ injected at the molecule-inlet port into each frozen
field (the last 20 ms time-averaged).
- The three cross-section conditions from handoff section 3 (commit
  5bca552). The table below is the energy-dependent default.
- Both constant-sigma bounds give ~0% extracted at every flow from 4 SCCM
  up.
- Ions were added until >= 150 were extracted per flow, and >= 600 across
  15-21 SCCM (+-8% at 95%).

**Thermalized** means < 0.1 eV kinetic energy at the box exit. That cut sits
in an empty gap between two populations:
- Gas-entrained cold ions, 10-30 meV. The lab-frame energy is mostly the
  ~120-175 m/s forward drift in the neon jet.
- Ions that never cooled, 1-400 eV.

### Results (energy-dependent sigma)

| Flow (SCCM) | ions traced | extracted, all (95% CI) | extracted, thermalized | back out inlet | median pump-out time | median exit energy |
|---|---|---|---|---|---|---|
| 1 | 19000 | 8.55% (8.16-8.95) | 0.00% | 0% | 3.4 us | 178 eV |
| 2 | 19000 | 1.25% (1.10-1.42) | 0.04% | 0% | 9.1 us | 10 eV |
| 4 | 21000 | 0.72% (0.61-0.84) | 0.70% | 0% | 3.8 ms | 15 meV |
| 5 | 20000 | 1.02% (0.89-1.17) | 1.01% | 0% | 8.3 ms | 15 meV |
| 8 | 15000 | 1.03% (0.88-1.20) | 1.03% | 1% | 26.5 ms | 18 meV |
| 10 | 11000 | 1.45% (1.24-1.69) | 1.45% | 1% | 32.7 ms | 19 meV |
| 15 | 24600 | 2.60% (2.41-2.80) | 2.60% | 6% | 48.0 ms | 22 meV |
| 16 | 18800 | 3.21% (2.96-3.47) | 3.21% | 7% | 47.5 ms | 22 meV |
| **17** | 18800 | **3.48% (3.23-3.76)** | 3.48% | 8% | 46.7 ms | 22 meV |
| **18** | 18800 | **3.43% (3.17-3.70)** | 3.43% | 12% | 51.0 ms | 23 meV |
| 19 | 18800 | 3.30% (3.05-3.56) | 3.30% | 16% | 51.0 ms | 23 meV |
| 20 | 19200 | 3.20% (2.96-3.46) | 3.20% | 21% | 50.2 ms | 23 meV |
| 21 | 26700 | 2.45% (2.27-2.64) | 2.45% | 25% | 55.5 ms | 24 meV |
| 22 | 6300 | 2.65% (2.28-3.08) | 2.65% | 27% | 58.8 ms | 23 meV |
| 23 | 5900 | 2.83% (2.44-3.29) | 2.83% | 35% | 50.4 ms | 25 meV |
| 25 | 8700 | 2.09% (1.81-2.41) | 2.09% | 42% | 54.0 ms | 25 meV |
| 30 | 30850 | 0.53% (0.46-0.62) | 0.53% | 72% | 70.3 ms | 25 meV |
| 58 | 10200 | 0.00% (0.00-0.04) | 0.00% | 99% | - | - |

### Reading it

- **1-2 SCCM:** the "extraction" is hot ions coasting straight through in
  microseconds, not a beam.
- **From 4 SCCM up:** every extracted ion leaves thermalized, carried out by
  the neon flow over milliseconds.
- **Peak:** extraction rises from ~1% (4-8 SCCM) to a broad maximum of
  3.2-3.5% at 16-20 SCCM, highest at 17-18 SCCM (statistically tied).
- **Past the peak:** extraction falls steeply. The outflowing neon increasingly
  blows the ions back out the molecule-inlet port: 8% at 17 SCCM, 72% at 30,
  99% at 58. So the optimum is set by that competition, not by diffusion to
  the walls.
- **21 SCCM** sits ~2 sigma below 20 and 22. It is likely a statistical or
  frozen-field fluctuation rather than structure.
- **Pump-out time** (injection to outlet) grows from ~4 ms at 4 SCCM to
  ~50 ms at the peak.

**Figures** (not in git; results/ is ignored):
- results/figures/fluor3d_scan_{extraction,pumpout,energy}.png and
  _summary.csv, from `tools/plot_fluor3d_scan.py`.
- A trajectory GIF of 40 ions at 17 SCCM in results/figures/fluor3d_traj/,
  from `tools/render_fluor3d_trajectories.py`.

**Pipeline** (in tools/wsl/):
- sbatch_fluor3d_scan.sh: neon run, self-chaining, preemption-safe.
- sbatch_fluor3d_scan_post.sh: asymmetry check + field freeze.
- sbatch_fluor3d_tracer_chunk.sh.
- Login-node drivers fluor3d_{scan,stats,peak}_driver.sh.

### Computed cross-section table (2026-09-27, commit 55b5f6b)

**Setup.**
- Same frozen neon fields, re-traced with the new tracer default
  `--sigma-model table`: the computed BaF+-Ne sigma_t(E_cm) in
  tracer/data/bafplus_ne_sigma_tracer.tsv.
- Three variants: central (the main condition), low and high. Low and high
  bracket the short-range wall, which matters above ~1 eV.
- Flows: 5, 10, 15 and 20 SCCM as requested, plus central at 25 and 30 SCCM
  to locate the peak.
- Statistics: >= 300 extracted ions per flow for central, >= 150 for low
  and high.
- `--sigma-model edep` under the new code reproduces the earlier scan
  bit-for-bit, so the table above still stands for that model.

| Flow (SCCM) | table central (95% CI; extracted/traced) | table low | table high | edep (previous model) | central: back out inlet | central: median pump-out |
|---|---|---|---|---|---|---|
| 5 | 3.20% (2.9-3.6; 320/10000) | 5.89% | 2.11% | 1.02% | 0% | hot ions, us |
| 10 | 4.58% (4.2-5.0; 550/12000) | 6.90% | 3.02% | 1.45% | 0% | 17.7 ms |
| 15 | 6.48% (5.8-7.2; 311/4800) | 8.67% | 4.25% | 2.60% | 1% | 39.2 ms |
| 20 | **8.62% (7.9-9.5; 414/4800)** | 12.58% | 6.65% | 3.20% | 5% | 43.9 ms |
| 25 | **8.56% (7.7-9.5; 308/3600)** | - | - | 2.09% | 8% | 50.6 ms |
| 30 | 5.60% (4.9-6.4; 224/4000) | - | - | 0.53% | 24% | 65.6 ms |

**What changes relative to edep:**
- **Extraction is 2-4x higher at every flow.** The central table peaks at
  ~8.6% over 20-25 SCCM, vs edep's ~3.5% at 17-18 SCCM, so the peak moves
  up in flow.
- **Mechanism: stopping depth.** sigma_t is much smaller than edep at high
  energy (4.9e-20 m^2 at 500 eV lab, vs edep's 1.5e-19 floor), so injected
  ions stop deeper in the chamber: ~14 vs 8.8 mm in the commit's synthetic
  check.
  - Deeper, they are farther from the molecule-inlet port, so fewer are
    blown back out: 5% at 20 SCCM (table) vs 21% (edep).
  - They are also closer to the outlet flow.
  - The low variant (weakest short-range wall) penetrates deepest and
    extracts the most (12.6% at 20 SCCM). The high variant extracts least
    (6.7%).
- **The thermal end is the same in both models** (2.0e-18 m^2), so from
  10 SCCM up essentially all extracted ions leave thermalized in every
  variant.
- **At 5 SCCM only 35% of central-table "extracted" ions are thermalized**
  (1.1% of injected). The rest coast straight through while still hot: the
  smaller high-energy sigma lets more ions stay fast at low density.
  Thermalized-only fractions at 5 SCCM: low 0.50%, central 1.12%,
  high 1.74%.
- **The spread is set by the short-range wall:** central 8.6%, bracketed by
  6.7-12.6% at 20 SCCM. It is not set by statistics. No BaF+-Ne cross
  section has been measured; Ba+ stands in for BaF+.

**Figures:**
- results/figures/fluor3d_scan_table_{extraction,pumpout,energy}.png, from
  `python tools/plot_fluor3d_scan.py --set table`.
- The extraction plot shows the three table variants, the central
  thermalized-only curve, and edep as a gray dashed reference.

### Thermal-injection test (2026-09-27)

**Setup.** BaF+ is injected already at the neon temperature: a 23 K
Maxwellian, no 500 eV beam (`--vx 0 -T 23`).
- Cross section: a single constant sigma = 1.93e-18 m^2.
- Neon fields: the same frozen 5-30 SCCM fields.
- Start points: three depths on the port axis
  (sbatch_fluor3d_tracer_chunk.sh COND=thermal-const1.93e-18[-x42|-mid]):
  - Inlet channel exit (x = 31.8mm). The literal "comes in thermalized".
    The outer port face at 30.48mm would send half the ions straight back
    out.
  - 10mm into the chamber (x = 41.8mm), about where 500 eV ions stop.
  - Chamber centre (x = 54.5mm), like molecules born inside the cell.

| Flow (SCCM) | centre: extracted (95% CI) | centre: median pump-out | 10mm in: extracted | 10mm in: back out inlet | inlet exit: extracted | inlet exit: back out inlet | 500 eV, table central | 500 eV, edep |
|---|---|---|---|---|---|---|---|---|
| 5 | 1.01% (0.90-1.13) | 16 ms | 0.25% | 4% | 0% (<0.04) | 95% | 3.20% | 1.02% |
| 10 | 5.43% (4.90-6.02) | 25 ms | 1.66% | 3% | 0% (<0.06) | 100% | 4.58% | 1.45% |
| 15 | 12.62% (11.4-14.0) | 30 ms | 3.58% | 4% | 0% (<0.05) | 100% | 6.48% | 2.60% |
| 20 | 22.96% (21.4-24.6) | 29 ms | 8.28% | 4% | 0% (<0.05) | 100% | 8.62% | 3.20% |
| 25 | 26.83% (25.1-28.6) | 32 ms | 10.42% | 2% | 0% (<0.05) | 100% | 8.56% | 2.09% |
| 30 | **34.30% (32.3-36.4)** | 37 ms | 10.65% | 5% | 0% (<0.06) | 100% | 5.60% | 0.53% |

**Reading it:**
- **Starting at the inlet, nothing gets out.**
  - The neon leaving through the molecule-inlet port sweeps a thermal ion
    straight back out, within ~0.1 ms.
  - A 500 eV ion only escapes this because it coasts several mm into the
    chamber before cooling.
- **So the injection depth, not the cross section, dominates extraction.**
  - A thermal ion 10mm in behaves much like the 500 eV edep/table runs:
    about 8-10% at 20-30 SCCM.
  - A thermal ion at the chamber centre is extracted far more often: 23% at
    20 SCCM and still rising to 34% at 30 SCCM, with no peak in this range.
  - It starts away from the inlet backflow and close to the outlet flow.
- **Why the 500 eV curves peak:** the peak at ~17-25 SCCM comes from ions
  that stop near the inlet being blown back out at high flow. Ions
  starting deeper don't see that loss, so their extraction keeps climbing
  with flow as diffusion to the walls slows.
- **Pump-out times are shorter than the 500 eV runs** (16-37 ms, vs
  ~45-65 ms), because the ions start closer to the outlet.

**Figures:** results/figures/fluor3d_scan_thermal_{extraction,pumpout,energy}.png,
from `python tools/plot_fluor3d_scan.py --set thermal`.

### Thermal injection deeper in, and at higher flow (2026-09-29 overnight)

**Setup.** Same thermal condition as above: a 23 K Maxwellian with constant
sigma = 1.93e-18 m^2, on the port axis. Two changes:
- Two more start points. The chamber spans x = 31.75 mm (inlet channel end)
  to 77.22 mm (outlet aperture plane: a 3.175 mm hole on the same axis).
  - 3/4 of the way through: x = 65.9 mm (COND suffix `-x66`).
  - 2 mm upstream of the outlet aperture: x = 75.2 mm (`-x75`).
- Higher flows, to look for a peak:
  - New neon runs at 40, 45 and 50 SCCM, same recipe as the scan, 0.3 s each.
  - The existing 58 SCCM field.
  - The existing 127 SCCM field, as an upper bracket. It was dropped earlier
    only for the 500 eV runs.
- Statistics: at least 300 extracted per point where affordable, capped at
  10,000 ions traced.
  - 10mm-in at 127 SCCM stopped at 122 extracted: each ion takes minutes
    there.
  - Earlier points keep their counts from the 2026-09-27 run.

**Extracted % (95% Wilson interval); median pump-out time for extracted ions:**

| Flow (SCCM) | 2 mm before outlet | 3/4 through | centre | 10 mm in |
|---|---|---|---|---|
| 5 | 26.6% (24.7-28.6); 0.62 ms | 4.5% (4.1-5.0); 6.1 ms | 1.0% (0.9-1.1); 16 ms | 0.25% (0.2-0.3); 28 ms |
| 10 | 57.4% (54.6-60.2); 0.57 ms | 17.0% (15.6-18.6); 8.0 ms | 5.4% (4.9-6.0); 25 ms | 1.7% (1.5-1.8); 41 ms |
| 15 | 72.7% (68.6-76.5); 0.56 ms | 33.6% (30.7-36.7); 8.5 ms | 12.6% (11.4-14.0); 30 ms | 3.6% (3.2-4.0); 49 ms |
| 20 | 87.1% (83.8-89.8); 0.53 ms | 47.0% (43.8-50.1); 7.0 ms | 23.0% (21.4-24.6); 29 ms | 8.3% (7.5-9.2); 50 ms |
| 25 | 92.7% (90.0-94.7); 0.57 ms | 58.1% (55.0-61.2); 6.8 ms | 26.8% (25.1-28.6); 32 ms | 10.4% (9.5-11.5); 52 ms |
| 30 | 95.8% (93.3-97.3); 0.55 ms | 66.0% (62.6-69.2); 7.3 ms | 34.3% (32.3-36.4); 37 ms | 10.7% (9.7-11.6); 65 ms |
| 40 | 98.2% (97.2-98.9); 0.54 ms | 73.4% (70.6-76.0); 6.1 ms | 50.9% (47.4-54.3); 27 ms | 15.3% (13.9-16.7); 56 ms |
| 58 | 100% (99.0-100); 0.56 ms | 88.5% (85.0-91.3); 5.7 ms | 58.0% (54.0-61.9); 36 ms | 21.7% (19.6-23.9); 74 ms |
| 127 | 100% (99.0-100); 0.58 ms | 97.0% (94.8-98.3); 6.3 ms | 86.3% (83.4-88.9); 24 ms | 21.8% (18.6-25.4); 60 ms |

45 and 50 SCCM are pending their neon runs.

**Reading it:**
- **No extraction peak up to 127 SCCM for a thermal ion, from any start
  point.** Every curve rises with flow, or saturates:
  - Near the outlet: ~100% from 58 SCCM up.
  - 10mm in: flat at ~22% between 58 and 127 SCCM.
  - This matches the earlier reading. The 500 eV peak near 17-25 SCCM comes
    from ions stopping near the inlet and being blown back out. A thermal
    ion that starts downstream of that backflow only gains from more flow:
    diffusion to the walls slows, and the outflow toward the outlet gets
    stronger.
- **Start point matters more than flow.** At 20 SCCM, extraction runs 87%,
  47%, 23% and 8% going from 2mm before the outlet back to 10mm in.
  - Wall losses dominate everywhere except near the outlet.
  - Only the 10mm-in start loses ions back out of the inlet: 2-5% at
    5-30 SCCM, 6-7% at 40-58 and 21% at 127. That is what flattens its
    curve at the top.
- **Pump-out time is set by distance, not flow.**
  - 2mm before the outlet: ~0.55 ms at every flow.
  - 3/4 through: 6-8 ms.
  - Centre: 25-37 ms.
  - 10mm in: 50-75 ms.
- **Caveat.** Everything at or above ~10 SCCM is on the 0.41mm chamber grid,
  which is coarser than the neon mean free path there (shaded in the plot).
  At 23 SCCM and up it is well under-resolved, and 127 SCCM is only a
  qualitative bracket. The port jets are approximate, but the stagnant
  chamber state is fine.

**Pipeline notes.**
- sbatch_fluor3d_tracer_chunk.sh gained the `-x66` and `-x75` start points.
- fluor3d_stats_driver.sh gained env options:
  - TOPT: several ';'-separated option sets, rotated per chunk.
  - CHUNKSCALE, PERFLOW, MINCHUNK.
- fluor3d_scan_driver.sh: PREP_ONLY=1 stops at the frozen field.
- sbatch_fluor3d_scan.sh: EXCLUDE is passed through the whole chain. Use a
  comma-free range such as `node[2000-2027]`, because --export splits on
  commas.
- Scheduling lessons:
  - Preemptable/normal node20xx nodes ran a 60000-step neon window in
    >30 min, vs 6-11 min on node14xx/16xx/19xx. The neon chains exclude them.
  - Short 30-min jobs on preemptable and normal started within minutes,
    even though `sbatch --test-only` estimated 1-2 days.

# Original (leaky-field) results

Follow-up to `docs/hpc-fluor-cell-3d-molecules-handoff.md`, working through its
section 2 smoke test and section 4 validation checklist against the completed,
settled neon run (`cases/fluor-cell-3d`, MDOT=6.002164e-08 kg/s, FILLN=5.4e21,
FNUM=1e11, settled from ~step 290000 onward, frozen at step 360000).

## Bugs found and fixed in `tracer/ParticleTracing3D.jl`

The handoff doc's own testing-status note said to treat the first real run
against this repo's actual geometry/field as the real first end-to-end test.
It was -- the section 2 smoke test (`-n 20`) surfaced two real bugs, neither
caught by the file's synthetic-cube local testing:

1. **Exit positions weren't clipped to the actual crossing point.**
   `freePath` samples free-flight steps up to 1000m as a safety cap (this
   case's box is ~0.2m). When a sampled step overshot the wall mesh or the
   simulation box, `getCollision`/`propagate` reported the raw, un-clipped
   endpoint as the particle's final position instead of where it actually
   crossed. In the first smoke test this put roughly half of 20 particles
   1-289 meters from a 0.2m box -- a silent bad value, not a crash, so
   easy to miss without checking individual rows (which the doc's own
   validation checklist calls for, and which is what caught this).

   Fixed by having `segment_hits_triangle` return the intersection
   fraction `t` (it already computed this internally, via the standard
   Moller-Trumbore parametrization, and simply discarded it), threading
   that through `getCollision` (now returns `(code, t)` instead of a bare
   code), and clipping `xnext` to the true crossing point in `propagate`
   before returning. The clip is nudged a hair (1e-9, sub-angstrom on this
   scale) past the exact boundary so a later re-classification call on the
   same output position reliably re-detects it as outside/on the mesh,
   rather than landing exactly on a boundary where a plain `>`/`<`
   comparison can go either way.

2. **A dummy call broke once (1)'s interface changed.** `SimulateParticles`
   makes one throwaway call to `propagate` purely to measure the output
   tuple's length for array allocation, using an inline placeholder
   `getCollision` of `(x,y)->true` (a stand-in for "the very first
   step always immediately collides", not real physics). That still
   returned a bare `Bool` after the interface changed to `(code, t)`,
   causing a `BoundsError: attempt to access Bool at index [2]` the
   moment the real fix went in. Fixed to `(x,y)->(1,1.0)`, matching the
   new interface.

Both call sites at `getCollision(SPAWN_REF, xpart)` (the `--spawnclip`
check) and the exit-classification re-check were updated for the new
`(code, t)` return shape as well.

Verified after the fix: all exit positions land within the domain; the
[300-particle statistics runs](#results) below ran cleanly with no
further position anomalies.

## Validation checklist (handoff doc section 4)

| Check | Result |
|---|---|
| Geometry read counts | 1127 points / 2258 triangles -- exact match to `check_surf3d.py`'s own count for the same file |
| Reproducibility | Identical output across repeated runs at the same `--seed`/`--threads` |
| Energy sanity | Exit speed retains 81-94% of the 24843.71 m/s injection speed even after collisions -- not down to the neon field's own ~1-200 m/s scale. Not a bug: BaF+ (156u) on Ne (20u) is a heavy-on-light collision (max ~40% KE transfer for a dead-on hit, less for glancing), and mean collision count here is only 0.6-2.1 -- nowhere near enough hits for real thermalization in this short a transit. |
| Collision-count sanity | Path-integrated mean-free-path estimate along the inlet-to-outlet centerline, using the real settled neon field's own density, gives ~5.4 expected collisions at sigma=5.0e-19; observed mean was 2.06. Same order of magnitude; the gap is plausibly explained by particles scattering off the direct centerline after their first hit, shortening their remaining path through the dense gas region relative to the naive straight-line estimate. |
| sigma sensitivity | Run at both 5.0e-19 and 1.5e-19 m^2 -- the transmission fraction moves meaningfully between them (see below), so this cross-section uncertainty is not negligible for a precise number. |

**Classification methodology note**: the handoff doc's literal "transmitted
if `xnext` near 0.08255m [the outlet]" criterion does not hold in practice.
The outlet is an open channel, not a wall -- a particle that passes it with
positive vx just keeps coasting through the near-vacuum downstream region to
the simulation box's far edge (x=0.2m) rather than stopping near the outlet
itself. Classification used instead: a particle that exits at a box face
with x beyond the outlet counts as transmitted; exits at a box face with
x at or before the outlet count as lost back toward the inlet; anything
that does NOT end at a box face (i.e. an interior position) is a wall hit
(stuck). This is well-justified given the mesh is confirmed watertight: a
particle can only leave the housing through a port opening or a wall hit,
so a box-exit downstream of the outlet essentially had to have passed
through it.

## Results

n=300 per condition, `--seed 42`, traced against the one completed neon
run (flow rate corresponding to MDOT=6.002164e-08 kg/s / ~4 SCCM):

| sigma (m^2) | Transmitted | Lost (back out inlet) | Stuck on wall | Mean collides | Ballistic (0 collides) |
|---|---|---|---|---|---|
| 5.0e-19 | 72.3% | 0.0% | 27.7% | 2.06 | 16.0% |
| 1.5e-19 | 89.3% | 0.0% | 10.7% | 0.61 | 57.0% |

The transmission fraction moves by 17 percentage points across the
cross-section uncertainty band -- both values say "mostly transmitted"
qualitatively, but a precise number should be reported as a range, not
collapsed to one point estimate, given this spread.

Notably, "lost back out the inlet" was 0% at both cross sections --- no
particle backscattered hard enough to exit back out the inlet side of the
box at this flow rate. Only one flow rate (the completed neon run) was
traced; a flow-rate scan would need additional neon runs at different
MDOT/FILLN.
