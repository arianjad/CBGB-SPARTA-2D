# Handoff: fluor-cell-3d work, moving from ORCD Engaging to MIT subMIT

Written 2026-10-01 on ORCD by the Claude session that did this work. You are the
next Claude, on subMIT. This file is self-contained.

Read in this order:
1. This file.
2. Your memory notes. They were copied from ORCD; section 7 lists which parts are
   now stale.
3. `docs/fluor-cell-3d-molecules-findings.md`, only when you need the full
   history. It has every result table and the reasoning behind it.

`docs/MIGRATION-subMIT.md` holds the setup commands the user ran. Section 1 here
checks that they worked.

---------------------------------------------------------------------------

## 0. What the user expects from you

- Continue the physics work (section 5). First confirm the environment works
  (section 1), and update your memory notes for subMIT (section 7).
- Working style the user asked for (also in memory):
  - **Never `git push` without asking.** Local commits are fine.
  - During long or unattended runs, check every ~30 min with a
    `CronCreate` cron at off-minutes, e.g. `13,43 * * * *`. End each check with a
    **one-line status update**, and keep check-in replies short.
  - Watch disk quota (`quota -s` on subMIT). Delete only clearly useless files of
    this project, never other projects' files, and ask before deleting raw data.
  - When the user asks "why" about the physics, answer with measured numbers
    from the runs, not just theory. For example, they once asked why
    higher-energy ions didn't go deeper. The answer was measured stopping depths
    plus the ln E collision argument.
  - Run heavy work (parsing GB-size files, builds, simulations) as batch jobs,
    never on the login node. On subMIT, interactive work must stay under 15 min.
  - When asked to run overnight, keep the cluster busy. If the main job finishes
    early, pick a sensible follow-up and say what you chose.

## 1. Verify the environment (do this first, report results in one short list)

Paths: `W=/work/submit/$USER`, `C=/ceph/submit/data/user/<u>/$USER`.
The repo is `$W/CBGB-SPARTA-2D`, branch `experiments`.

```bash
cd $W/CBGB-SPARTA-2D
git log --oneline -3              # should include this handoff commit and 0414bb1 (streaming freeze)
ls -l ~/orcd/scratch              # symlink -> $C ; all scripts use ~/orcd/scratch/cbgb-sparta-3d
ls ~/orcd/scratch/cbgb-sparta-3d  # fluor-scan  fluor-scan-fine  fluor-scan-ports  fluor-cell-3d-ne-sealed
echo $SBATCH_PARTITION            # submit   (overrides the scripts' "#SBATCH --partition=mit_*")
ls $HOME/opt/sparta-27Aug2026/bin/spa_mpi      # SPARTA, built with SPARTA_BIGBIG=1
.local/julia-1.9.4/bin/julia --project=tracer -e 'using CSV, DataFrames; println("julia ok")'
source .venv/bin/activate && python -c "import numpy, pandas, pyvista; print('py ok')"
bash -c 'source tools/wsl/site_env.sh; echo $SITE'   # must print: submit
```

If anything is missing, fix it with the matching step in `docs/MIGRATION-subMIT.md`.

Then run two smoke tests as batch jobs. Use `sbatch -p submit`, keep them small,
and check `seff <jobid>` afterwards to learn subMIT's speed and memory.

1. **Tracer.** Re-trace 300 ions on the coarse 17 SCCM field:
   ```bash
   sbatch -p submit -t 00:30:00 --cpus-per-task=16 --mem=48G \
     --export=ALL,SCCM=17,SEED=9001,NPART=300,COND=edep tools/wsl/sbatch_fluor3d_tracer_chunk.sh
   ```
   - Output: `~/orcd/scratch/cbgb-sparta-3d/fluor-scan/sccm17/post/tracer_edep_s9001.out`.
   - Classify it with the snippet in section 3.4.
   - Expect ~3.5% extracted. With 300 ions the CI is wide; just check it isn't 0
     or absurd.
   - If Julia dies with "Illegal instruction" (rc 132), the precompile cache
     doesn't fit that CPU type: wipe `~/.julia/compiled` and precompile again.
     On ORCD this happened on some node types.
2. **Neon / SPARTA.** A 2000-step run of the coarse case in a throwaway directory:
   ```bash
   mkdir -p $C/smoke && cp cases/fluor-cell-3d/{in.ne_fluor3d_sealed_run,cell_fluor3d_sealed.surf,ne.species,ne.vss} $C/smoke/
   sbatch -p submit -n 48 -t 00:30:00 --mem=64G --wrap "cd $C/smoke && module load mpi && \
     mpirun -np 48 $HOME/opt/sparta-27Aug2026/bin/spa_mpi -in in.ne_fluor3d_sealed_run \
     -var MDOT 2.550920e-07 -var FILLN 2.2950e+22 -var FNUM 8.5e+11 -var DT 2.9e-07 \
     -var DUMPFREQ 1000 -var CHECKPOINT_STEPS 1000 -var NCHUNKS 2 > run.log 2>&1"
   ```
   - Expect ~1.49M grid cells, ~1M particles, and "Particles stuck = 0".
   - ORCD reference speed: ~630 s per 60000 steps on 40 cores (~0.01 s/step).
   - Delete `$C/smoke` afterwards.

## 2. Project in brief

- **What is modeled:** a 3D DSMC (SPARTA 27Aug2026) model of 23 K neon buffer gas
  in a sealed fluorescence/cooling cell, `cases/fluor-cell-3d/cell_fluor3d_sealed.surf`.
  - Coordinates in m.
  - The molecule inlet port face is at x = 30.48 mm and the chamber starts at
    31.75 mm.
  - The outlet aperture (3.175 mm) is at x = 77.22 mm, and its channel ends at
    82.55 mm.
  - All ports sit on the axis y = z = 19.05 mm.
  - The buffer gas enters at (34.29, 25.4, 0) mm.
- **The tracer:** BaF+ ions are flown through a time-averaged ("frozen") neon
  field by the Julia tracer `tracer/ParticleTracing3D.jl`.
- **Figure of merit:** the fraction extracted through the outlet, especially the
  thermalized fraction (< 0.1 eV at exit). Losses are wall hits and "back out
  the inlet".
- **Cross sections** (`--sigma-model`):
  - `table` (default): the computed BaF+-Ne table, `tracer/data/bafplus_ne_sigma_tracer.tsv`.
    It was extended to ~19 keV lab in 41e1cc2 and matches ZBL nuclear stopping
    to 1-2%. Variants: `--sigma-table-variant central|low|high`.
  - `edep`: the older 1/v model, floored at 1.5e-19 m^2. The user calls it "the
    1/v cross section".
  - `const`.
- **Flow scale:** flows are in SCCM. Chamber density scales ~linearly (aperture
  law), e.g. 2.15e22 m^-3 at 17 SCCM. The neon mean free path is ~0.4 mm x
  4/SCCM.

## 3. How to run things (the pipeline)

All scripts live in `tools/wsl/`, are submitted from the repo root
(`$SLURM_SUBMIT_DIR`), and keep their data under
`~/orcd/scratch/cbgb-sparta-3d/<RUNBASE>/sccm<F>/`, which is the /ceph symlink.
Each flow directory has its run files, `field.grid` (the raw dumps; large), and
`post/`. `post/` holds `frozen.grid`, `asymmetry.txt` and `tracer_<COND>_s<SEED>.out`.

### 3.1 Neon run: `sbatch_fluor3d_scan.sh`

- Restartable and self-chaining. Each job runs MAXCHUNKS averaging windows,
  checkpoints, then resubmits itself until TOTAL steps.
- A per-flow lock (`lock.d`) lets only one job work on a flow.
- Per-flow DT / window / TOTAL come from `fluor3d_flow_params.sh`.
- Options (all env, passed through the chain):

| Option | Meaning |
|---|---|
| `SCCM` | flow (required) |
| `MAXCHUNKS` | windows per job; 0 = run to the end in one job |
| `CHAIN_TIME` | walltime of resubmitted jobs |
| `GRIDREF` | chamber refinement per axis (8 = 0.41 mm, 16 = 0.20 mm); FNUM is scaled by (8/n)^3 |
| `RUNBASE` | data subdir (default `fluor-scan`) |
| `SIM_S` | simulated seconds (default 0.3) |
| `WINDOW` | averaging window = checkpoint interval, in steps |
| `INPUT=ports` | port-refined input (`cases/fluor-cell-3d/in.ne_fluor3d_ports_run/_resume`) |
| `PORTREF` | level-3 refinement at the ports |
| `FNUM_SET`, `DT_SET` | overrides |
| `CHAIN_NODES` | nodes per job |
| `EXCLUDE` | node exclude for the chain; pass a file path, not a comma list |

- NP (the task count) must stay fixed for the life of a run, because restart
  files are written per processor.
- **On subMIT** jobs can run up to 6 days, so prefer one long job, e.g.:
  ```bash
  sbatch -p submit -n 48 -t 2-00:00:00 --mem=128G \
    --export=ALL,SCCM=17,RUNBASE=<newbase>,MAXCHUNKS=0,CHAIN_TIME=2-00:00:00,<grid options> tools/wsl/sbatch_fluor3d_scan.sh
  ```
  The ORCD 15/30-min chaining was only to fit ORCD queues. Ask the user what
  job size and length they prefer on subMIT; fair-share there favors smaller
  requests.

### 3.2 Post: `sbatch_fluor3d_scan_post.sh`

- Writes `post/frozen.grid` (the last 20 ms averaged) with
  `tools/field3d_freeze_stream.py`.
- Run with `NPART=0` (field prep only).
- Pass the same `RUNBASE`, `WINDOW`, `DT_SET` as the neon run.
- `ASYM_FRAMES=0` skips the full asymmetry check and writes only the chamber
  summary to `asymmetry.txt`. Use it on big grids: `tools/fluor3d_asymmetry.py`
  loads all its frames at once and needed >400 GB on the 0.20 mm grid.
- The streaming freeze needs ~8.5 GB on 2.6M cells and ~1.5 GB/frame on 11.5M
  cells (4.5 min). Example:
  ```bash
  sbatch -p submit -t 01:00:00 --cpus-per-task=4 --mem=64G \
    --export=ALL,SCCM=17,NPART=0,RUNBASE=<base>,WINDOW=<w>,ASYM_FRAMES=0 tools/wsl/sbatch_fluor3d_scan_post.sh
  ```

### 3.3 Tracer: `sbatch_fluor3d_tracer_chunk.sh` + `fluor3d_stats_driver.sh`

- One chunk:
  ```bash
  sbatch -p submit ... --export=ALL,SCCM=<F>,SEED=<unique>,NPART=<n>,COND=<cond>[,RUNBASE=<base>] tools/wsl/sbatch_fluor3d_tracer_chunk.sh
  ```
  - COND: `edep`, `table-central|low|high`, `const2e-18`, `const1.5e-19`, or
    `thermal-const1.93e-18[-x42|-mid|-x66|-x75]` (a thermal start at depth).
  - Beam energy suffix: `-e<eV>`, e.g. `table-central-e5000`.
  - Seeds must be unique per COND per flow.
  - Output: `post/tracer_<COND>_s<SEED>.out`.
- The driver keeps topping up chunks until TARGET ions are extracted or NMAX
  traced. Its state (in-flight job ids, next seed, chunk size) lives in `post/`.
  Start it on the login node, detached:
  ```bash
  RUNBASE=<base> COND=<cond> FLOWS="17 20" TARGET=600 NMAX=40000 MAXJOBS=6 PERFLOW=6 CHUNKSCALE=0.3 \
  TOPT="-p submit -t 01:00:00 --cpus-per-task=16 --mem=96G" \
  setsid nohup bash tools/wsl/fluor3d_stats_driver.sh >> driver_<name>.log 2>&1 < /dev/null &
  ```
  - ALWAYS pass TOPT on subMIT; the built-in default names an ORCD partition.
  - TOPT may hold several `;`-separated option sets, which are rotated per chunk.
  - A timed-out chunk halves the chunk size (down to MINCHUNK).
- Check progress via `post/.stats_jobs*` and the tracer outputs, not only the
  log: after restarts the log can contain NUL padding that hides lines.

### 3.4 Analysis and plots

- Classify tracer outputs:
  ```python
  import glob, sys; sys.path.insert(0, "tools")
  from tracer_classify3d import load, classify
  from plot_fluor3d_scan import wilson
  r = load(glob.glob(P + "/tracer_edep_s*.out")); tr, lost, wall = classify(r)
  p, lo, hi = wilson(tr.sum(), len(r))
  ```
  - Columns: idx x y z xn yn zn vx vy vz collides time.
  - Kinetic energy at exit for extracted ions: `0.5*156.325*1.66054e-27*|v|^2`.
    Thermalized means < 0.1 eV.
  - Use numpy argsort, not `DataFrame.sort_values`, which raises IndexError in
    the Python 3.14/pandas venv.
- Plots:
  - `python tools/plot_fluor3d_scan.py --set {edep,table,thermal,energy,kev}`
    writes `results/figures/fluor3d_scan_<set>_*`.
  - `python tools/plot_fluor3d_energy20.py`.
  - Trajectory GIFs: `tools/render_fluor3d_trajectories.py` with
    `--sccm --energy --linear-after`. Its input is a tracer run with
    `--saveall 1 --trajprint 40` (see git log 3b600bd/15d049c for the job
    pattern). pyvista renders off-screen (`PYVISTA_OFF_SCREEN=true`). Check it
    works on subMIT compute nodes.
- Before writing plotting or figure code, load the `dataviz` skill. The reference
  palette slots are used throughout.

## 4. Results so far (details and CIs are in the findings doc)

- **500 eV injection, edep:** extraction peaks ~17-20 SCCM at 3.2-3.5%. At higher
  flow, inlet backflow wins.
- **Thermal-start test** (23 K ions, constant 1.93e-18):
  - Extraction depends mainly on where the ion starts. At 20 SCCM: 8% (10 mm in),
    23% (chamber centre), 47% (3/4 through), 87% (2 mm before the outlet).
  - No peak up to 127 SCCM. Denser gas slows diffusion to the walls.
- **Stopping depth vs energy:**
  - With edep, ions stop ~7 mm in for 0.5-1 keV (stopping ~ ln E, ~20% energy
    lost per collision).
  - With the table, at 20 SCCM: 12 / 20 / 26 / 35 mm at 0.5 / 2 / 5 / 10 keV.
- **keV injection with the table:**
  - **Best: 5 keV at ~40 SCCM, 26% thermalized** (22-26% over 25-50 SCCM).
  - 2 keV peaks ~15% (20-30 SCCM); 500 eV table peaks 8.6%.
  - At <= 10 SCCM keV ions fly out hot.
  - At 20 SCCM / 5 keV, the sigma low/high variants give 13.6 / 16.9% vs 19.4%
    central.
- **Caveats:** not modeled: electronic stopping (LSS ~30-50% of nuclear at 1-10
  keV), BaF+ dissociation in keV collisions, space charge. The grid is
  under-resolved at >= 10 SCCM, and section 5 tests that.

## 5. Current work: grid convergence at 17 SCCM (500 eV, edep)

| Grid | RUNBASE | Cells | Chamber n (m^-3) | Extracted (95% CI) | Back out inlet |
|---|---|---|---|---|---|
| COARSE: 0.41 mm chamber | fluor-scan | 1.5M | 2.152e22 | 3.48% (3.23-3.76), 655/18800 | 8.1% |
| PORTS: 0.41 mm + ~0.135 mm boxes at the inlet, gas jet and outlet (3 levels, volume-weighted, PORTREF=3, FNUM_SET=1.2052e22 per unit volume, DT 2e-7) | fluor-scan-ports | 2.6M | 1.915e22 | 3.74% (3.47-4.04), 633/16920 | 2.1% |
| FINE: uniform 0.20 mm chamber (GRIDREF=16, DT 2.9e-7, WINDOW 8000, 0.15 s) | fluor-scan-fine | 11.5M | 1.985e22 | tracer was running on ORCD at handoff (see below) | (pending) |

What it shows so far:
- Extraction is insensitive to resolution within the CIs.
- Inlet blow-back and chamber density are sensitive. Finer grids lower the
  chamber density by 8-11%: the outlet conductance is better resolved. This is
  not a weighting artifact; FINE uses no weighting and agrees.
- Finer ports cut blow-back from 8.1% to 2.1%. Those ions end on walls instead.

FINE status:
- If `fluor-scan-fine/sccm17/post/tracer_edep_s*.out` has >= 600 extracted ions
  in the copied data, the result is complete.
  - The ORCD session was meant to write it up in the findings doc ("Grid
    convergence at 17 SCCM") and rebuild the package. Check `git log`.
  - If that doc section is missing, write it.
- If FINE outputs are absent or partial, the frozen field
  (`post/frozen.grid`, 1.5 GB) is there if the final package was copied.
  - Re-trace with the driver (section 3.3, `RUNBASE=fluor-scan-fine`,
    `TOPT` --mem >= 96G) to 600 extracted.
  - If even frozen.grid is missing, the raw dumps stayed on ORCD. Ask the user.

Next steps (propose them, ask before big runs):
1. Finish FINE and the doc section, then commit locally.
2. **0.10 mm port grid** (`INPUT=ports PORTREF=4`). It needs SPARTA built with
   `SPARTA_BIGBIG=1`; on ORCD it failed with "Hierarchical grid induces cell IDs
   that exceed 32 bits".
   - FNUM_SET stays per unit volume (1.2052e22 at 17 SCCM).
   - Use DT ~1.5e-7.
   - Smoke-test first. Check the buffer-gas emission rate equals
     mflow*dt/(m*fnum*V_cell) per step (it did at PREF=3).
3. Re-check the best case (5 keV, ~40 SCCM, table) on a refined grid, since
   backflow is grid-sensitive.
4. Sanity checks the user agreed were worthwhile (not done yet):
   - Tracer diffusion coefficient in uniform static gas vs Chapman-Enskog and
     measured Ba+ mobility in Ne.
   - Flow-following in uniform flowing gas.
   - Range vs SRIM/LSS including electronic stopping.
   - Chunk-to-chunk binomial variance check.
   - Neon mass balance.

## 6. Lessons learned (still apply)

- Never edit a driver `.sh` while an instance of it is running. Bash reads
  scripts incrementally, and running drivers died mid-loop.
- Kill by PID. Never `pkill -f`: it killed the calling shell.
- A wait loop must not `pgrep -f` a pattern that appears in its own command line.
- `--export` splits values on commas. Pass node lists as a file path
  (`--exclude=<file>`) or a comma-free range.
- Check new sbatch options with `sbatch --test-only`, but its start estimates
  were often pessimistic.
- SPARTA volume weighting (`global weight cell volume`) makes FNUM per unit cell
  volume. Exterior 3.25 mm cells then get very few particles, so the exterior
  plume is noisy.
- Raw dumps are big: 1.55-2.8 GB per frame at 11.5M cells. Resumed jobs re-dump
  an empty duplicate frame at their start step; the freeze tools skip those.
- Low-flow keV ions leave hot. Report thermalized yields.
- Julia "Illegal instruction" (rc 132) on some CPU types means a stale precompile
  cache.
- Detached processes (`setsid nohup`) survived some session restarts on ORCD but
  not all. After any restart, check drivers are alive and re-launch them with
  the same env. Their state is in the post dirs, so nothing is lost.

## 7. Memory notes to update on subMIT

These ORCD facts are now wrong or irrelevant:
- partitions `mit_quicktest/mit_preemptable/mit_normal`, account `mit_general`,
  and the QOS caps (quicktest 48 CPUs);
- bad-node lists (node17xx/18xx/21xx, node20xx/25xx, node1617);
- the 15/30-min job chaining;
- `~/orcd/.quota`;
- the `/etc/claude-code/CLAUDE.md` ORCD policy.

Add the subMIT facts:
- partition `submit` (~1000 CPUs, up to 6 days), `submit-gpu`;
- `module load mpi`;
- quotas: home 10 GB, /work 100 GB, /ceph 1 TB; check with `quota -s`;
- code lives in /work, data in /ceph behind the `~/orcd/scratch` symlink;
- interactive work < 15 min;
- fair-share favors small requests.

Keep the physics results and the user preferences.
</content>
</invoke>
