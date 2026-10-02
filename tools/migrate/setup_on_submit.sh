#!/bin/bash
# Unpack the ORCD migration package on subMIT and put everything in place.
#   bash setup_on_submit.sh [TARGET]      (default TARGET=/work/submit/$USER/buffer_gas)
# Run it from the package directory (where CBGB-SPARTA-2D.bundle is).
# It is safe to re-run: existing pieces are left alone.
#   - clones the repo INTO TARGET (TARGET becomes the repo root; start Claude there)
#   - unpacks results/figures, the Julia tarball, the post-processed data
#   - makes ~/orcd/scratch -> /ceph/submit/data/user/<u>/$USER (all scripts use that path)
#   - copies Claude's memory notes to ~/.claude/projects/<TARGET path with non-alphanumerics -> '-'>/memory
# Builds (Julia packages, Python venv, SPARTA) are listed at the end; run them as batch jobs.
set -euo pipefail
PKG="$(cd "$(dirname "$0")" && pwd)"
TARGET="${1:-/work/submit/$USER/buffer_gas}"
C="/ceph/submit/data/user/${USER:0:1}/$USER"
[ -f "$PKG/CBGB-SPARTA-2D.bundle" ] || { echo "run this from the package directory" >&2; exit 1; }

echo "== repo -> $TARGET"
mkdir -p "$TARGET"
if [ -d "$TARGET/.git" ]; then echo "repo already there"
elif [ -z "$(ls -A "$TARGET")" ]; then
  git clone -q "$PKG/CBGB-SPARTA-2D.bundle" "$TARGET"
  git -C "$TARGET" checkout -q experiments
  git -C "$TARGET" remote set-url origin git@github.com:heericmc/sparta.git
else echo "ERROR: $TARGET is not empty and is not a git repo" >&2; exit 2; fi
git -C "$TARGET" log --oneline -1

echo "== figures, Julia"
[ -d "$TARGET/results/figures" ] || tar -C "$TARGET" -xf "$PKG/results-figures.tar"
mkdir -p "$TARGET/.local"
[ -d "$TARGET/.local/julia-1.9.4" ] || tar -C "$TARGET/.local" -xzf "$PKG/julia-1.9.4-linux-x86_64.tar.gz"
cp -n "$PKG/julia-1.9.4-linux-x86_64.tar.gz" "$TARGET/.local/" || true

echo "== data -> $C/cbgb-sparta-3d, and ~/orcd/scratch -> $C"
mkdir -p "$C/cbgb-sparta-3d"
[ -d "$C/cbgb-sparta-3d/fluor-scan" ] || tar -C "$C/cbgb-sparta-3d" -xf "$PKG/cbgb-sparta-3d-post.tar"
mkdir -p "$HOME/orcd"
[ -e "$HOME/orcd/scratch" ] || ln -s "$C" "$HOME/orcd/scratch"
ls -ld "$HOME/orcd/scratch"

echo "== Claude memory + settings"
KEY=$(echo "$TARGET" | sed 's/[^A-Za-z0-9]/-/g')
K="$HOME/.claude/projects/$KEY"
mkdir -p "$K/memory"
cp -rn "$PKG/claude/memory/." "$K/memory/"
[ -f "$HOME/.claude/settings.json" ] || cp "$PKG/claude/settings.json" "$HOME/.claude/settings.json"
ls "$K/memory"

echo "== Slurm partition for the scripts' '#SBATCH --partition=mit_*' lines"
grep -q "SBATCH_PARTITION=submit" "$HOME/.bashrc" || echo 'export SBATCH_PARTITION=submit' >> "$HOME/.bashrc"

cat <<EOF

Done. Files are in place. Remaining builds (each as a batch job, from $TARGET):
  cd $TARGET
  # Julia tracer packages
  sbatch -p submit -c 4 --mem=16G -t 01:00:00 --wrap ".local/julia-1.9.4/bin/julia --project=tracer -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'"
  # Python env (conda python 3.14 + venv; install Miniforge into /work/submit/$USER/miniforge3 first if you have no conda)
  conda create -y -p /work/submit/$USER/conda/cbgb python=3.14 && /work/submit/$USER/conda/cbgb/bin/python -m venv .venv \\
    && source .venv/bin/activate && pip install -r requirements-3d-freeze.txt
  # SPARTA with 64-bit cell IDs, installed where the scripts expect it
  sbatch -p submit -c 8 --mem=16G -t 01:00:00 --wrap "module load mpi; SPARTA_BIGBIG=1 RUN_TESTS=0 BUILD_JOBS=8 \\
    SPARTA_SRC_ROOT=/work/submit/$USER/src/sparta-27Aug2026 SPARTA_BUILD_ROOT=/work/submit/$USER/build/sparta-27Aug2026 \\
    bash tools/wsl/build_sparta_plain.sh"
Then: source ~/.bashrc; cd $TARGET; claude   -- and paste the first message from docs/MIGRATION-subMIT.md.
EOF
