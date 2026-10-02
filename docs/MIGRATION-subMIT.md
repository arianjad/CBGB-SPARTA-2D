# Moving the fluor-cell-3d work from ORCD Engaging to MIT subMIT

Built on ORCD by `tools/migrate/make_migration_package.sh` into
`~/orcd/scratch/cbgb-migrate/` (this file is that package's README.md).

Target layout on subMIT (user heeric):
- repo root: `/work/submit/heeric/buffer_gas` (start Claude here)
- Claude memory: `~/.claude/projects/-work-submit-heeric-buffer-gas/memory/`
- data: `/ceph/submit/data/user/h/heeric/cbgb-sparta-3d`, reached through the
  symlink `~/orcd/scratch` so all scripts work unchanged

## Package contents (~6 GB)

| File | Contents |
|---|---|
| `CBGB-SPARTA-2D.bundle` | the whole git repo: all branches/commits, including commits never pushed to GitHub |
| `cbgb-sparta-3d-post.tar` | every flow's `post/`: frozen neon fields, all BaF+ tracer outputs, asymmetry checks; run inputs and scan logs. NOT the raw `field.grid` dumps or restarts (~280 GB, they stay on ORCD) |
| `results-figures.tar` | `results/figures` + `results/exports` |
| `julia-1.9.4-linux-x86_64.tar.gz` | Julia for the tracer |
| `claude/` | memory notes, `settings.json`, old session transcripts (`transcripts.tgz`) |
| `setup_on_submit.sh` | puts everything above in place on subMIT |

## Step 1: copy the package to subMIT

Option A runs on subMIT and needs no new keys. You'll get a password/Duo prompt
for ORCD:
```bash
mkdir -p /ceph/submit/data/user/h/heeric/cbgb-migrate
rsync -av --partial --progress heeric@orcd-login.mit.edu:orcd/scratch/cbgb-migrate/ \
      /ceph/submit/data/user/h/heeric/cbgb-migrate/
```
If the ORCD path doesn't resolve, use the absolute path
`/orcd/scratch/orcd/013/heeric/cbgb-migrate/`. Rerun the same command to resume
an interrupted copy.

Option B is a push from ORCD. First upload the ORCD public key
`~/.ssh/id_ed25519_submit.pub` at https://submit-portal.mit.edu/, then on ORCD:
```bash
rsync -av --partial --progress -e "ssh -i ~/.ssh/id_ed25519_submit" \
  ~/orcd/scratch/cbgb-migrate/ heeric@submit.mit.edu:/ceph/submit/data/user/h/heeric/cbgb-migrate/
```

## Step 2: unpack and place everything (on subMIT)

```bash
cd /ceph/submit/data/user/h/heeric/cbgb-migrate
bash setup_on_submit.sh /work/submit/heeric/buffer_gas
```

The script:
- clones the repo into `/work/submit/heeric/buffer_gas`, branch `experiments`;
- unpacks figures, Julia and the data;
- creates the `~/orcd/scratch` link;
- copies the memory notes to `~/.claude/projects/-work-submit-heeric-buffer-gas/memory/`;
- adds `export SBATCH_PARTITION=submit` to `~/.bashrc`.

It's safe to re-run. At the end it prints the three build commands below.

## Step 3: builds (from /work/submit/heeric/buffer_gas)

```bash
cd /work/submit/heeric/buffer_gas && source ~/.bashrc
# Julia tracer packages
sbatch -p submit -c 4 --mem=16G -t 01:00:00 --wrap ".local/julia-1.9.4/bin/julia --project=tracer -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'"
# Python env: conda python 3.14 + venv. subMIT recommends conda under /work;
# install Miniforge into /work/submit/heeric/miniforge3 first if you have no conda.
conda create -y -p /work/submit/heeric/conda/cbgb python=3.14
/work/submit/heeric/conda/cbgb/bin/python -m venv .venv
source .venv/bin/activate && pip install -r requirements-3d-freeze.txt
# SPARTA with 64-bit cell IDs (allows the 0.10 mm port grid), installed at
# $HOME/opt/sparta-27Aug2026 where the scripts expect it
sbatch -p submit -c 8 --mem=16G -t 01:00:00 --wrap "module load mpi; SPARTA_BIGBIG=1 RUN_TESTS=0 BUILD_JOBS=8 \
  SPARTA_SRC_ROOT=/work/submit/heeric/src/sparta-27Aug2026 SPARTA_BUILD_ROOT=/work/submit/heeric/build/sparta-27Aug2026 \
  bash tools/wsl/build_sparta_plain.sh"
```

If those builds aren't done yet, the new Claude can run them too. They are in
its checklist.

## Step 4: Claude Code

```bash
curl -fsSL https://claude.ai/install.sh | bash   # if not installed; then /login inside claude
cd /work/submit/heeric/buffer_gas && claude
```

Paste this as the first message:

```text
We just moved this project from the ORCD Engaging cluster to MIT subMIT. The previous
Claude session (on ORCD) wrote a detailed handoff for you.

1. Read docs/HANDOFF-subMIT.md completely before doing anything else. It is your
   runbook: what the project is, how to run every pipeline stage on subMIT, the
   results so far, the work in progress, lessons learned, and my preferences.
   Your memory notes were copied over from ORCD too. Some of their cluster
   details are now wrong (the handoff's section 7 lists which).

2. Verify the setup with the checklist in handoff section 1: repo/branch, the
   ~/orcd/scratch -> /ceph symlink, SBATCH_PARTITION=submit, SPARTA (built with
   64-bit cell IDs), Julia + tracer packages, the Python .venv, and site_env.sh
   detecting "submit". Finish any build in docs/MIGRATION-subMIT.md step 3 that
   isn't done. Then run the two small smoke tests (tracer and neon) as batch
   jobs. Report the results to me as a short list.

3. Update your memory notes for subMIT: replace the ORCD-specific facts
   (partitions, account, bad-node lists, ~/orcd/.quota, short job chaining) with
   the subMIT ones from the handoff. Keep the physics results and my
   preferences.

4. Then continue the grid-convergence work (handoff section 5). Check whether the
   uniform 0.20 mm (FINE) tracer result came over complete. If not, finish it
   from the frozen field. Write the "Grid convergence at 17 SCCM" doc section if
   it isn't there yet, and commit locally. After that, propose the next step:
   I'm most interested in the 0.10 mm port-refined grid now that SPARTA can be
   built with 64-bit IDs. Ask me before launching big runs, and ask how large
   and long I want subMIT jobs to be.

Ground rules: never git push without asking; run heavy work only as batch jobs
(interactive work under 15 minutes on subMIT); watch disk quota (quota -s);
during long runs check every 30 minutes and give me one-line updates.
```

## What stays on ORCD

The raw neon `field.grid` dumps (coarse flow scan ~160 GB, uniform 0.20 mm run
~180 GB) and the Takahashi validation run, in
`~/orcd/scratch/cbgb-sparta-3d/`. They are only needed to re-freeze over a
different time window or to make field plots. ORCD scratch is not backed up.
</content>
</invoke>
<invoke name="Bash">
<parameter name="command">[ -f ~/.ssh/id_ed25519_submit ] || ssh-keygen -q -t ed25519 -N "" -C "heeric@orcd-for-submit" -f ~/.ssh/id_ed25519_submit; cat ~/.ssh/id_ed25519_submit.pub; cd ~/CBGB-SPARTA-2D && git add docs/MIGRATION-subMIT.md docs/HANDOFF-subMIT.md tools/migrate/ && git commit -q -m "Migration: setup_on_submit.sh, target /work/submit/heeric/buffer_gas

setup_on_submit.sh clones the bundle into the target folder (the repo root
where Claude starts), unpacks figures/Julia/data, links ~/orcd/scratch to
/ceph, and copies the memory notes to the matching
~/.claude/projects/<path-key>/memory. README/handoff updated for that
layout, with pull (from subMIT) and push (from ORCD) copy options and the
full first message for the new session.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" && git log --oneline -1 | cat && bash tools/migrate/make_migration_package.sh 2>&1 | grep -E "setup_on|README|ready"