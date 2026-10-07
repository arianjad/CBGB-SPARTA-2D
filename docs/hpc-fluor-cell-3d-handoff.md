# Handoff: running the 3D fluorescence cooling cell on HPC

From-scratch guide, no prior context on the repo assumed. Goes from a clean
checkout to a completed buffer-gas-only 3D neon DSMC run.

Scope: buffer gas (neon) only. Nothing involving the BaF+ ion/molecule tracer.
That is a separate, larger job (extending `tracer/ParticleTracing.jl` from
2D-axisymmetric to true 3D) that starts only after the buffer-gas results below
have been reviewed and the flow conditions confirmed. Stop and report back
after this doc.

## 0. What this case is

`cases/fluor-cell-3d/` is a true 3D SPARTA case (unlike the axisymmetric
`cases/b5-lean/` and `cases/takahashi-single/`) for a real CAD geometry: a
cryogenic buffer-gas beam cell with a side-mounted gas inlet, so the flow is
non-axisymmetric and can't be reduced to 2D. The surface mesh came from a STEP
CAD file, defeatured (mounting screw holes plugged, leaving the 3 functional
ports), and checked against the source CAD for volume, watertightness and
normal orientation (details in the commit that added this case).

The three ports, as SPARTA surf types in `cell_fluor3d.surf`:

| Type | What | Notes |
|---|---|---|
| 1 | Generic wall | 23K diffuse-reflecting, everywhere else |
| 2 | Buffer gas inlet | Sealed: a capped stub, not an open hole. Gas is injected via `emit/surf` mass-flow, same convention as `cases/b5-lean` and `cases/takahashi-3d`. |
| 3 | Molecule inlet channel | 2.54mm, open to exterior. SPARTA doesn't track molecules; tagged only so its neon flux can be measured separately. |
| 4 | Molecule outlet channel | 3.175mm, open to exterior. Same. |

Buffer gas is 23K neon, not the helium used by the other cases.
`ne.species`/`ne.vss` are new (no neon data in the repo before);
`ne.vss` uses literature VHS parameters (Bird 1994, Table A1: d=2.79e-10m,
omega=0.66), not the simplified hard-sphere convention `he.vss` uses.

## 1. Environment setup

Follow the main [README](../README.md) section "1. Install and clone" for
the base environment (WSL2/Ubuntu if on Windows, apt dependencies, Julia;
Julia isn't needed for this pass but installing it now saves a step later).
Use this repo's URL, not the upstream `arianjad/CBGB-SPARTA-2D` one:

```bash
mkdir -p ~/code
cd ~/code
git clone https://github.com/heericmc/sparta.git
cd sparta
git checkout experiments
git pull

python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-2d.txt
```

Build the pinned SPARTA solver as the README describes. This case needs
nothing beyond the standard plain-MPI build:

```bash
bash tools/wsl/build_sparta_plain.sh --print-plan
BUILD_JOBS=4 TEST_JOBS=1 bash tools/wsl/build_sparta_plain.sh
export SPARTA_PLAIN_EXE="$HOME/opt/sparta-27Aug2026/bin/spa_mpi"
```

If the pinned build (`27Aug2026`) is already installed on this machine, just set `SPARTA_PLAIN_EXE` to it and skip the build.

## 2. Sanity-check the surf file itself (no SPARTA needed)

Confirm the geometry file transferred correctly:

```bash
cd cases/fluor-cell-3d
python3 check_surf3d.py
```

Expect: `3469 points, 6842 triangles`, triangle counts by type
`{1: 6684, 2: 30, 3: 64, 4: 64}`, `Free edges: 0`, and `OK -- watertight,
ready for read_surf` at the end. If any of those differ, stop: the file didn't
transfer correctly (check line endings / binary-vs-text transfer mode) or was
edited.

## 3. Smoke test: does SPARTA accept the geometry?

A `run 0`: only checks that SPARTA reads the surface without errors (format,
watertight from SPARTA's side, consistent normals). Takes seconds.

```bash
cd ~/code/sparta/cases/fluor-cell-3d
"$SPARTA_PLAIN_EXE" -in in.smoketest_3d_fluor
```

Expect it to run to completion and print stats for a single step with no
errors. If `read_surf` errors out (e.g. about unclosed surfaces or bad
normals), stop and report the full error. Don't edit the surf file by hand to
get around it.

## 4. The real buffer-gas run

`in.ne_fluor3d_scoping` is the actual scoping deck, structured like
`cases/takahashi-3d/in.he_takahashi_3d_scoping`. It needs 5 values supplied
on the command line via `-var` (not hardcoded, same as every other case here;
see `tools/run_helium.sh`):

| Variable | Meaning | Starting value (see derivation below) |
|---|---|---|
| `MDOT` | buffer gas mass flow, kg/s | `6.002164e-08` |
| `FILLN` | initial prefill density, m^-3 | `5.4e21` |
| `FNUM` | simulation particle weight | `2.5e17` (from `cases/b5-lean`, similar density regime; check particles per cell once running, see step 5) |
| `DT` | timestep, s | `1e-7` (from `cases/b5-lean` as a starting point, not verified for this grid, see step 5) |
| `STEPS` | total steps | start with a short check, `2000`, before a long run |

Where `MDOT` and `FILLN` come from (redo this if the target flow changes; both
scale with SCCM):

```python
import math
kB = 1.380649e-23  # J/K
amu = 1.66053906660e-27  # kg
m_ne = 20.1797 * amu
T = 23.0
v_mean_ne = math.sqrt(8*kB*T/(math.pi*m_ne))          # 155.34 m/s

sccm = 4.0                        # <-- change this to scan flow rate
ndot = sccm * 4.478e17            # particles/s per SCCM (standard conversion)
mdot = ndot * m_ne                # kg/s  -> this is MDOT

RAP = 3.175e-3 / 2                # molecule-outlet radius, the dominant exit
filln = 4*ndot / (1.07*v_mean_ne*math.pi*RAP**2)   # m^-3 -> this is FILLN
```

Same aperture-law method as `cases/b5-lean/gen_b5.py`'s He fill estimate,
redone for neon at 23K. Caveat for this geometry: the B5 cell has one exit
aperture, this cell has two open channels (molecule inlet, 2.54mm, and
molecule outlet, 3.175mm). The formula only counts the outlet, so it
overestimates; the real equilibrium density will be lower. It doesn't need to
be exact: `FILLN` only changes how long the run takes to settle, not the
equilibrium itself (set by the geometry, `MDOT` and boundary conditions).
Same caveat as the README's He case ("This is an initial estimate, not a
measured equilibrium density").

Run the short check first:

```bash
cd ~/code/sparta/cases/fluor-cell-3d
"$SPARTA_PLAIN_EXE" -in in.ne_fluor3d_scoping \
  -var MDOT 6.002164e-08 -var FILLN 5.4e21 -var FNUM 2.5e17 \
  -var DT 1e-7 -var STEPS 2000
```

If that completes cleanly, do a real settling run with much larger `STEPS`.
The He example needed 120,000 steps (12ms simulated) for 10ms settling + 2ms
sampling at its `DT`. This geometry and grid haven't been characterized, so
the same step count, `DT` and grid may not be right; see step 5 before trusting
a long run.

For multi-rank runs, put `mpirun -np N` in front, as for every other case:

```bash
mpirun -np 4 "$SPARTA_PLAIN_EXE" -in in.ne_fluor3d_scoping -var MDOT ... [etc]
```

## 5. Known gaps

This case is less mature than `cases/b5-lean`/`cases/takahashi-3d`:

- `create_box`/`create_grid` in `in.ne_fluor3d_scoping` are placeholders,
  not derived from a mean-free-path calculation. The mature 2D case
  (`in.he_b5_mflow`) uses nested 2:1 refinement regions sized so the finest
  cells resolve the smallest geometric features (see its `region level*`
  blocks and `gen_b5.py`'s `nonface()` checks, which assert cell boundaries
  don't land exactly on a geometric feature). This case has a single flat
  grid level. Before trusting quantitative results, check
  the actual mean free path in the settled field against the cell size (see
  [Controlling numerical accuracy](accuracy.md) for the general method) and
  refine the grid around the three ports (1.588/2.54/3.175mm, all much
  smaller than a flat grid cell) the same way the 2D case does.
- No 3D equivalent of `tools/check_grid_2d.py` (the 2D nested-refinement
  checker). Write one or check by hand.
- No 3D equivalent of `tools/plot_fields_b5.py` for the field dump. That tool
  only reads the 2D `Points`/`Lines` surf format (`tools/plot_fields.py`'s
  `wall_segments()` parses `Lines` records; this surf file has `Triangles`). For now, inspect `field.grid` as text, or write a
  minimal 3D reader before trying to reuse the existing plotting code.
- No run launcher like `tools/run_helium.sh`/`run_lean.sh` (those are
  2D-specific: they stage `.surf`/`.species` files and write a provenance
  manifest). This doc calls `spa_mpi` directly. A launcher is worth writing
  if there are many runs, but not needed for a first result.

None of these block a first run, but they matter before trusting a
quantitative result (a specific density or flow number).

## 6. What to check afterward, and what to send back

Once you have a completed run:

1. Confirm it finished: check the tail of the run's stdout/log for
   the final `run` command's stats, and that SPARTA exited without an error
   (no crash, no `ERROR:` lines).
2. Look at the `dump gd grid all ... field.grid` output: the neon
   density/velocity/temperature field, dumped every 20000 steps (averaged
   over 2000-step windows per `fix ag ave/grid`, same cadence as the He
   cases). Density near the sealed gas-inlet cap should be higher than
   downstream (gas flowing from the type-2 cap toward the two open
   channels), and the last couple of frames should agree (settled, not
   transient).
3. Report: whether it ran cleanly, steps needed to settle, the
   `MDOT`/`FILLN`/`FNUM`/`DT`/grid actually used (if changed from the
   starting values above), and either the raw
   `field.grid` file or a summary of the settled density/velocity near each
   of the three ports. That gets reviewed before choosing a flow rate and
   moving to the molecule/ion tracer work.

The BaF+/ion tracer extension waits for that review; this doc stops here.
