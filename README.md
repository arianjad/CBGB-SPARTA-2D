# CBGB-SPARTA-2D

Simulate helium flow in a cryogenic buffer-gas cell, then follow molecules
through that helium field. This repository contains the **2D axisymmetric**
SPARTA inputs, a Julia molecule tracer, and Python analysis and plotting tools.
It runs independently of the research repository it was extracted from.

The included B5 example is a 4 K cell with a 38.1 mm inner diameter,
63.5 mm body length, and a 5 mm diameter aperture through a 1.5 mm plate.
Helium enters at 4 SCCM. The example molecule source is hot BaF; its mass,
collision cross section, temperature, and birth distribution are configurable.

**Start with the five steps below.** They take you from installation to helium
maps, molecular extraction/speed/angle summaries, and trajectory figures.
See [Gotchas](docs/2d-gotchas.md) when something looks wrong.

![Twenty recorded BaF paths over the helium flow field](docs/images/trajectories.png)

Example output from the tutorial: these first 20 molecules all hit a wall.
The full 1,000-molecule run is used for the extraction statistics.

## 1. Install and clone

Use Linux, or an Ubuntu terminal in WSL2 on Windows. The reference setup is
Ubuntu 22.04 with Python 3.10 and Julia 1.9.4. Native Windows and macOS are
not covered by these shell scripts. On Windows, install WSL first from
[Microsoft's instructions](https://learn.microsoft.com/windows/wsl/install),
then run the following **inside Ubuntu**, not PowerShell.

```bash
sudo apt update
sudo apt install -y git curl ca-certificates build-essential cmake \
  openmpi-bin libopenmpi-dev python3 python3-venv python3-dev
mkdir -p ~/code
cd ~/code
git clone https://github.com/arianjad/CBGB-SPARTA-2D.git
cd CBGB-SPARTA-2D

python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-2d.txt

bash tools/install_julia.sh
export JULIA="$PWD/.local/julia-1.9.4/bin/julia"
JULIA_NUM_PRECOMPILE_TASKS=4 "$JULIA" --project=tracer -e 'using Pkg; Pkg.instantiate()'
```

The Julia installer supports Linux x86-64 and ARM64, checks the official
archive checksum, and installs only under this clone's ignored `.local/`.
If Julia 1.9.4 is already installed, skip the installer and set `JULIA` to its
executable. Keep `tracer/Manifest.toml`; it pins the tracer packages. The first
Julia invocation compiles packages and is slower than later runs.

Build the pinned SPARTA solver:

```bash
bash tools/wsl/build_sparta_plain.sh --print-plan
BUILD_JOBS=4 TEST_JOBS=1 bash tools/wsl/build_sparta_plain.sh
export SPARTA_PLAIN_EXE="$HOME/opt/sparta-27Aug2026/bin/spa_mpi"
```

This downloads SPARTA, compiles it, runs upstream tests, and installs without
sudo. The recipe pins `27Aug2026` / `95b9abaa8bd548991cc3c3f1c58b34722f7ade74`;
older builds lack the inlet's `mflow` option. If that exact build is already
installed, set `SPARTA_PLAIN_EXE` to it instead. Existing build/install
directories are protected; use fresh `SPARTA_BUILD_ROOT` and
`SPARTA_INSTALL_ROOT` values when rebuilding.

## 2. Set up the example

Run from the repository root. Repeat these exports in each new terminal:

```bash
source .venv/bin/activate
export JULIA="$PWD/.local/julia-1.9.4/bin/julia"
export SPARTA_PLAIN_EXE="$HOME/opt/sparta-27Aug2026/bin/spa_mpi"
export DSMC_RUN_ROOT="$PWD/results/he"
export DSMC_LOG_ROOT="$PWD/results/logs"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python cases/b5-lean/gen_b5.py
```

The generator checks the geometry/station coordinates and regenerates the
five tracked surface files. An unchanged example should leave `git diff`
empty. Run names below are examples: choose a **new name** for every attempt,
including after a failed run. Keep each completed run directory intact.

## 3. Compute the helium field

```bash
bash tools/run_helium.sh first-he
```

The default is four MPI ranks and **12 ms of simulated time**: 120,000 steps
at 100 ns per step. The cell, feed stub, and aperture channel start prefilled;
the downstream plume starts empty. The fill uses the analytic aperture-law
estimate printed by the geometry generator, `FILLN=2.3438e21 m^-3` at 4 SCCM.
This is an initial estimate, not a measured equilibrium density.

Allow the first **10 ms for settling**. The deck saves a 2 ms averaged field
every 20,000 steps; the final frame therefore samples 10–12 ms and is the
frame used for molecule tracing. Ten milliseconds is a conservative starting
prescription for this tutorial, not proof that every observable or new
geometry has equilibrated. If visible transients persist, extend the run
before interpreting results. Equilibrium investigations are outside this
getting-started example.

For a quick installation check before the full example:

```bash
STEPS=2000 RANKS=1 bash tools/run_helium.sh install-check
```

That check is only 0.2 ms and does **not** produce a usable molecule field.
Use a rank count your machine can support; don't run competing solver jobs.
The launcher refuses an existing run directory or a second active SPARTA job.
Wall-clock duration depends on the machine and particle count; 12 ms is
physical time, not a runtime estimate.

Check completion and inspect the grid:

```bash
cat "$DSMC_RUN_ROOT/first-he/rc.sentinel"  # must be 0
python -m json.tool "$DSMC_RUN_ROOT/first-he/manifest.json"
python tools/check_grid_2d.py "$DSMC_RUN_ROOT/first-he/field.grid"
```

The manifest must say `complete`. `run.log` holds the solver output; errors
or `failed` status should be resolved before launching molecules.

## 4. Trace molecules and make plots

```bash
N=1000 THREADS=4 SEED=42 TRAJPRINT=20 \
  bash tools/run_2d_standalone.sh \
  "$DSMC_RUN_ROOT/first-he" "$PWD/results/molecules/first-baf"
```

This converts the final helium field, traces 1,000 molecules, scores extraction
and crossings at the observation plane, and saves figures. `TRAJPRINT=20`
records paths for the first 20 molecules; use `0` (the default) to omit path
recording. Increase `N` for smaller molecular sampling error; recorded paths
are illustrative and not an unbiased sample of extracted molecules alone.

The default [source file](cases/b5-lean/baf-hot-source.conf) uses BaF
(`156.325 u`, He collision cross section `2.7e-18 m^2`) with a 1000 K isotropic
translational Maxwellian, born uniformly in a 12.7 mm radius sphere tangent
to the cell side wall at mid-length. Molecules stick when they hit a wall.

To try another molecule/source, copy that file and pass the copy as the third
argument. All lengths are metres; `SPAWN_SIZE_M` is a radius for `uniformball`
and a Gaussian standard deviation for `gaussball`. Use species-appropriate
cross sections; the BaF value is not a universal molecule–He constant.

```bash
cp cases/b5-lean/baf-hot-source.conf my-source.conf
# Edit my-source.conf, then use a new output name:
N=1000 THREADS=4 bash tools/run_2d_standalone.sh \
  "$DSMC_RUN_ROOT/first-he" "$PWD/results/molecules/my-source" my-source.conf
```

## 5. Read and keep your results

Open `results/molecules/first-baf/figures/` for helium maps, molecule phase
space, wall-hit/source maps, and (when requested) trajectory figures.

| Output in the molecule directory | What it contains |
|---|---|
| `exit.json`, `common.json` | Exit statistics and first-crossing acceptance at the configured plane |
| `trace/legs.out`, `trace/recs.csv`, `trace/spawn.csv` | Terminal legs, crossing records, and molecule birth positions |
| `field/` | The exact converted helium field, wall geometry, and field provenance |
| `source.conf`, `commands.sh` | A saved copy of the source and the executed commands |
| `manifest.json`, `rc.sentinel` | Input hashes, software identity, stage/status, and final return code |

Successful runs have sentinel `0` and manifest status `complete`. Keep **both**
helium and molecule directories: the molecule receipt refers back to the
helium run. Generated runs, environments, and downloaded tools are ignored by
Git. A snapshot's time is in `field/provenance.json` (`source_timestep=120000`
for this example).

`common.json` reports illustrative aperture/angle/time cuts. State the
observation plane, radius, and angular cut when quoting an acceptance.
Its binomial error describes finite molecule sampling on one frozen helium
field; it excludes helium-field noise and numerical/model uncertainty.

## Changing the simulation

The helium launcher accepts `MDOT` (kg/s), `FILLN` (m^-3), `FNUM`, `DT` (s),
`STEPS`, `RANKS`, and `SEED` through environment variables. Changing the flow
does not automatically change the fill. Recalculate it with the aperture-law
expression in `cases/b5-lean/gen_b5.py`; reassess particle count, grid, and
timestep as the physical conditions change. Smaller `FNUM` means more
simulated helium particles at fixed geometry and density.

The geometry generator and `in.he_b5_mflow` must be edited together for a new
cell: refinement regions and diagnostic stations are geometry-dependent.
See [Gotchas](docs/2d-gotchas.md) before changing them. The helium field is
axisymmetric; the molecule positions and velocities are Cartesian in 3D.
Tracer axial `z` corresponds to SPARTA `x`; tracer radial distance
`sqrt(x^2+y^2)` corresponds to SPARTA `y`.

## Tests, scope, and credits

Run `bash tools/test.sh` after setup for portable geometry, launcher, and
analysis checks. The short helium check plus the full molecule command above
exercise the installed solver/tracer path. See [validation](docs/validation.md)
for the release checks and their limits.

This is a one-way-coupled, frozen-helium-field transport model with elastic
collisions and a prescribed molecular source. It does not model rotational
cooling, chemistry, or ablation feedback on helium. The tracer samples local
fields at the start of each free flight and uses a mean-speed approximation
for the collision rate. Installation success alone does not validate a
quantitative molecular prediction.

The tracer derives from Cal Miller's
[DSMC_Simulations](https://github.com/cal-miller-harvard/DSMC_Simulations).
SPARTA comes from [sparta/sparta](https://github.com/sparta/sparta).
Original material uses the MIT license; the tracer and derived code retain
GPL-3 terms. See [third-party notices](THIRD_PARTY_NOTICES.md).
