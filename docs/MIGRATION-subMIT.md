# Moving CBGB-SPARTA-2D (fluor-cell-3d work) from ORCD Engaging to MIT subMIT

Prepared 2026-10-01 on ORCD. On ORCD, `tools/migrate/make_migration_package.sh`
builds a migration package in `~/orcd/scratch/cbgb-migrate/`. This file is that
package's README. Replace `<user>` with your subMIT username and `<u>` with
its first letter.

## What is in the package (~6 GB)

| File | Contents |
|---|---|
| `CBGB-SPARTA-2D.bundle` | The whole git repo, every branch and commit, including commits never pushed to GitHub |
| `cbgb-sparta-3d-post.tar` | Post-processed simulation data. Every flow's `post/` (frozen neon fields `frozen.grid`, all BaF+ tracer outputs, asymmetry checks), the 4 SCCM reference, and the grid-convergence runs. Also the run inputs and `scan.log`s. NOT included: the raw `field.grid` dumps (~250 GB) and restart files; they stay on ORCD scratch. |
| `results-figures.tar` | `results/figures` + `results/exports` (plots, summary CSVs, GIFs) |
| `julia-1.9.4-linux-x86_64.tar.gz` | The Julia version the tracer uses |
| `claude/` | Claude Code memory notes, settings.json, past session transcripts |

## Step 1: subMIT access (once)

1. Make sure you can `ssh <user>@submit.mit.edu`. If you can't, upload a public key at
   https://submit-portal.mit.edu/.
2. To copy straight from ORCD, the ORCD account needs a key subMIT accepts. On ORCD:
   `ssh-keygen -t ed25519 -f ~/.ssh/id_ed25519_submit`, then upload
   `~/.ssh/id_ed25519_submit.pub` on the subMIT portal too.

## Step 2: copy the package (run on ORCD, ~6 GB)

```bash
# on orcd-login (rsync resumes if interrupted: just rerun it)
rsync -av --partial --progress -e "ssh -i ~/.ssh/id_ed25519_submit" \
  ~/orcd/scratch/cbgb-migrate/ <user>@submit.mit.edu:/ceph/submit/data/user/<u>/<user>/cbgb-migrate/
```

## Step 3: set up subMIT (run on submit.mit.edu)

```bash
C=/ceph/submit/data/user/<u>/<user>      # 1 TB, data
W=/work/submit/<user>                    # 100 GB, code + software (home is only 10 GB)

# 3a. repo (keep code in /work, not /ceph)
cd $W && git clone $C/cbgb-migrate/CBGB-SPARTA-2D.bundle CBGB-SPARTA-2D
cd CBGB-SPARTA-2D && git checkout experiments
git remote set-url origin git@github.com:heericmc/sparta.git   # needs a GitHub key on subMIT (deploy key or your own)
tar -xf $C/cbgb-migrate/results-figures.tar                    # -> results/figures, results/exports

# 3b. data: every script uses ~/orcd/scratch/cbgb-sparta-3d, so point that at /ceph
mkdir -p $C/cbgb-sparta-3d && tar -C $C/cbgb-sparta-3d -xf $C/cbgb-migrate/cbgb-sparta-3d-post.tar
mkdir -p ~/orcd && ln -s $C ~/orcd/scratch

# 3c. Slurm: the scripts' "#SBATCH --partition=mit_*" lines are overridden by this
echo 'export SBATCH_PARTITION=submit' >> ~/.bashrc && source ~/.bashrc

# 3d. Julia (the tracer)
mkdir -p .local && tar -C .local -xzf $C/cbgb-migrate/julia-1.9.4-linux-x86_64.tar.gz
.local/julia-1.9.4/bin/julia --project=tracer -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# 3e. Python (post-processing/plots); subMIT recommends conda under /work
#     (install Miniforge into $W/miniforge3 first if you have no conda)
conda create -y -p $W/conda/cbgb python=3.14
$W/conda/cbgb/bin/python -m venv .venv
source .venv/bin/activate && pip install -r requirements-3d-freeze.txt

# 3f. SPARTA, with 64-bit cell IDs (allows the finer 0.10 mm port grid);
#     build as a batch job, not on the login node. The scripts expect
#     $HOME/opt/sparta-27Aug2026/bin/spa_mpi (7 MB, fine in home).
sbatch -p submit -c 8 --mem=16G -t 01:00:00 --wrap "module load mpi; \
  SPARTA_BIGBIG=1 RUN_TESTS=0 BUILD_JOBS=8 SPARTA_SRC_ROOT=$W/src/sparta-27Aug2026 \
  SPARTA_BUILD_ROOT=$W/build/sparta-27Aug2026 bash tools/wsl/build_sparta_plain.sh"
```

The batch scripts detect subMIT on their own (`tools/wsl/site_env.sh`: `module load mpi`
there, the ORCD modules here). Node-exclude lists in the instructions name ORCD nodes;
drop them on subMIT.

## Step 4: Claude Code on subMIT

```bash
curl -fsSL https://claude.ai/install.sh | bash      # then run `claude` and /login
mkdir -p ~/.claude && cp $C/cbgb-migrate/claude/settings.json ~/.claude/
# Memory is stored per working directory: the folder name is that directory's path
# with "/" replaced by "-". Start Claude from the repo, so:
K=~/.claude/projects/$(echo $W/CBGB-SPARTA-2D | tr / -)
mkdir -p $K && cp -r $C/cbgb-migrate/claude/memory $K/
# optional, old transcripts (only needed to scroll back through old sessions):
# tar -C $K -xzf $C/cbgb-migrate/claude/transcripts.tgz
```

Then `cd $W/CBGB-SPARTA-2D && claude`, and start with:

> We just moved this project from ORCD Engaging to subMIT. Read
> docs/HANDOFF-subMIT.md and your memory notes, check the setup, update the
> memory notes for subMIT, then continue from where the handoff leaves off.

## What stays on ORCD

- Raw neon `field.grid` dumps (coarse flow scan ~160 GB, uniform 0.20 mm grid-convergence
  run ~100 GB) and the Takahashi validation run: `~/orcd/scratch/cbgb-sparta-3d/`.
  Only needed to re-freeze over a different time window or to make field plots.
- ORCD scratch is not backed up. Copy anything you want to keep before leaving.
</content>
</invoke>
