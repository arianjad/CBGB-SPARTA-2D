# Handoff: fluor-cell-3d work, moving from ORCD Engaging to MIT subMIT (2026-10-01)

You are picking this up on subMIT after the user moved the project. First read
`docs/MIGRATION-subMIT.md` (setup steps), then this file, then the memory notes
(copied from ORCD; they describe ORCD and need updating). The full narrative and
every result table are in `docs/fluor-cell-3d-molecules-findings.md`.

## The project in one paragraph

This is a 3D DSMC (SPARTA) model of neon buffer gas flowing through a sealed
fluorescence/cooling cell (`cases/fluor-cell-3d`). Gas enters through a buffer-gas
inlet, BaF+ ions enter through a molecule inlet port, and gas leaves through an
outlet aperture, all on one axis.
- A Julia tracer (`tracer/ParticleTracing3D.jl`) flies BaF+ ions through a
  time-averaged ("frozen") neon field.
- The figure of merit is the fraction of ions extracted through the outlet,
  especially the thermalized fraction (< 0.1 eV at exit).
- Pipeline:
  - `tools/wsl/sbatch_fluor3d_scan.sh`: neon run per flow, a self-chaining
    restartable batch job.
  - `sbatch_fluor3d_scan_post.sh`: freezes the last 20 ms into
    `post/frozen.grid` and writes `post/asymmetry.txt`.
  - `sbatch_fluor3d_tracer_chunk.sh`: one batch of tracer ions.
  - `fluor3d_stats_driver.sh`: login-node driver that tops up tracer chunks until a
    target number of extracted ions.
  - `tools/plot_fluor3d_scan.py --set {edep,table,thermal,energy,kev}`: plots.
- Cross sections:
  - `table`: the computed BaF+-Ne table, the default. Extended to ~19 keV lab in
    41e1cc2; it matches ZBL nuclear stopping to 1-2%.
  - `edep`: the older 1/v model, floored at 1.5e-19 m^2.

## Key results so far

- **500 eV injection, edep model:** extraction peaks ~17-20 SCCM at 3.2-3.5%,
  limited by inlet backflow at higher flow.
- **Thermal-start test:** extraction depends mainly on where ions stop. At 20 SCCM:
  8% from 10 mm in, 23% from the chamber centre, 47% from 3/4 through, 87% from
  2 mm before the outlet. No flow peak up to 127 SCCM.
- **Energy:**
  - With edep, ions stop ~7 mm in at any energy (stopping ~ ln E).
  - With the table, depth grows ~sqrt(E): 12 / 20 / 26 / 35 mm at 0.5 / 2 / 5 /
    10 keV (20 SCCM).
- **Best case so far: 5 keV at ~40 SCCM, 26% thermalized** (flat 22-26% over
  25-50 SCCM). 2 keV peaks ~15%. At <= 10 SCCM keV ions fly out hot.
  - The low/high sigma variants at 20 SCCM / 5 keV give 13.6 / 16.9% vs 19.4%
    central.
- **Caveats:** not modeled: electronic stopping (LSS ~30-50% of nuclear at 1-10
  keV), BaF+ dissociation, space charge. The grid is under-resolved at >= 10
  SCCM.

## Grid-convergence work in progress (the user's current focus)

The user wants better grid resolution. The test case is 17 SCCM, 500 eV, edep.

| Grid | Cells | Chamber n (m^-3) | Extracted (95% CI) | Back out inlet |
|---|---|---|---|---|
| COARSE: 0.41 mm chamber (original) | 1.5M | 2.152e22 | 3.48% (3.23-3.76) | 8.1% |
| PORTS: 0.41 mm + ~0.135 mm boxes at inlet/gas-jet/outlet (3 levels, volume-weighted) | 2.6M | 1.915e22 | 3.74% (3.47-4.04) | 2.1% |
| FINE: uniform 0.20 mm chamber | 11.5M | see below | see below | see below |

- **Runs:** `fluor-scan-ports/sccm17`, `fluor-scan-fine/sccm17` (same `post/` layout
  as the scan).
- **Inputs:**
  - `cases/fluor-cell-3d/in.ne_fluor3d_ports_run/_resume`.
  - `REFINE` in `in.ne_fluor3d_sealed_run`.
  - The scan-script env options GRIDREF/RUNBASE/SIM_S/WINDOW/INPUT/PORTREF/FNUM_SET/DT_SET
    are documented at the top of `sbatch_fluor3d_scan.sh` (commit 13674de).
- **PORTS:**
  - Extraction is unchanged within CI.
  - Finer ports cut inlet blow-back sharply (8.1% -> 2.1%), but those ions are
    lost to walls instead.
  - Chamber density is 11% lower. FINE also trended ~5% lower than COARSE, with
    no weighting, so this looks like a real resolution effect (better-resolved
    outlet conductance), not a weighting bug.
- **FINE status at handoff:**
  - The neon run finished 0.15 s (520000 steps).
  - Its post job OOM'd at 180 GB; a retry at 700 GB was running on ORCD.
  - The tracer had not started when the package was made.
  - If the user re-ran `tools/migrate/make_migration_package.sh` after FINE
    finished, `fluor-scan-fine/sccm17/post/` will hold `frozen.grid` +
    `tracer_edep_s*.out`. Otherwise its result is still on ORCD, or redo it on
    subMIT (raw dumps were not copied).
- **To finish:**
  - FINE's row of the table (density from `post/asymmetry.txt`, extraction via
    `tools/tracer_classify3d.py`).
  - A doc section "Grid convergence at 17 SCCM" in the findings doc.
  - A local commit.
- **Natural next steps** (ask the user first):
  - A 0.10 mm port grid (`PORTREF=4`) now that SPARTA can be built with 64-bit
    cell IDs (`SPARTA_BIGBIG=1`). On ORCD it failed with "Hierarchical grid
    induces cell IDs that exceed 32 bits".
  - Then re-check the 5 keV / 40 SCCM best case on a refined grid.

## First things to do on subMIT

1. Check the setup from MIGRATION-subMIT.md:
   - `git log -1` shows the handoff commit.
   - `~/orcd/scratch` -> /ceph link exists.
   - `$HOME/opt/sparta-27Aug2026/bin/spa_mpi` exists (built with BIGBIG).
   - `.local/julia-1.9.4/bin/julia --project=tracer` loads.
   - `.venv` has numpy/pyvista.
2. Smoke tests, each as a short batch job (subMIT: interactive work < 15 min,
   no heavy work on login nodes):
   - Re-trace ~200 ions with COND=edep on
     `~/orcd/scratch/cbgb-sparta-3d/fluor-scan/sccm17/post/frozen.grid`.
     Expect ~3.5% extracted.
   - A 2000-step neon run (the smoke-test pattern is in the findings doc / git
     history).
3. Update the memory notes. They hold ORCD facts that are now wrong:
   - partitions mit_quicktest/preemptable/normal, account mit_general;
   - bad-node lists (node17xx/18xx/21xx, node20xx/25xx);
   - the 15/30-min job chaining;
   - `~/orcd/.quota`.
   subMIT facts:
   - partition `submit`, up to 6 days per job, fair-share favors smaller jobs;
   - quotas: home 10 GB, /work 100 GB, /ceph 1 TB (`quota -s`);
   - no `.quota` file;
   - `module load mpi` for OpenMPI.
   With 6-day limits the neon chains can run in one long job: pass
   `MAXCHUNKS=0 CHAIN_TIME=<long>`. Ask the user about job sizes first.
4. The stats driver's built-in TOPT default still names mit_preemptable.
   On subMIT always pass `TOPT="-p submit -t <time> --cpus-per-task=16 --mem=96G"`,
   or edit the default. It is safe to edit while no driver is running.

## User preferences (also in memory)

- Never push to GitHub without asking. Local commits are fine.
- During unattended runs, check every 30 min and give 1-line status updates.
- Short walltimes were preferred on ORCD for queue priority. Ask how they want
  it on subMIT.
- Don't fill the disk; check quota. Delete only clearly useless project files,
  never other projects.
- When asked "why" about physics, answer with measured numbers from the runs.
</content>
</invoke>
