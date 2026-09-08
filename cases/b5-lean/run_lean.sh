#!/usr/bin/env bash
# run_lean.sh RUNSUB NPROC DECKNAME VARSTRING BINARY CASEDIR FILE...
# Portable lean runner. Invoke this tracked file from a Linux/WSL clone.
# Override roots/aliases with DSMC_REPO_ROOT, DSMC_RUN_ROOT, DSMC_LOG_ROOT,
# SPARTA_PLAIN_EXE.
set -u

if [ "$#" -lt 7 ]; then
  echo "usage: $0 RUNSUB NPROC DECKNAME VARSTRING BINARY CASEDIR FILE..." >&2
  exit 2
fi

sub=$1; np=$2; deck=$3; vars=$4; bin=$5; casedir=$6; shift 6
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
repo=${DSMC_REPO_ROOT:-$(cd "$script_dir/../.." && pwd -P)}
run_root=${DSMC_RUN_ROOT:-$HOME/runs/he}
log_root=${DSMC_LOG_ROOT:-$HOME/logs}
python_cmd=${DSMC_PYTHON:-python3}
manifest_tool=${DSMC_MANIFEST_TOOL:-$repo/tools/lean_manifest.py}
plain_exe=${SPARTA_PLAIN_EXE:-$HOME/opt/sparta-27Aug2026/bin/spa_mpi}
sparta_tag=${SPARTA_TAG:-}
sparta_commit=${SPARTA_COMMIT:-}
solver_receipt=${SPARTA_BUILD_RECEIPT:-}
runner_thread_budget=${DSMC_THREAD_BUDGET:-16}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
export NUMEXPR_NUM_THREADS=${NUMEXPR_NUM_THREADS:-1}

case "$sub" in ""|*/*|.|..) echo "REFUSE: RUNSUB must be one directory name" >&2; exit 2 ;; esac
case "$casedir" in ""|/*|../*|*/../*|*/..) echo "REFUSE: CASEDIR must stay under cases/" >&2; exit 2 ;; esac
case "$np" in ''|*[!0-9]*|0) echo "REFUSE: NPROC must be a positive integer" >&2; exit 2 ;; esac
case "$runner_thread_budget" in ''|*[!0-9]*|0) echo "REFUSE: DSMC_THREAD_BUDGET must be a positive integer" >&2; exit 2 ;; esac
per_rank_threads=1
for count in "$OMP_NUM_THREADS" "$OPENBLAS_NUM_THREADS" "$MKL_NUM_THREADS" "$NUMEXPR_NUM_THREADS"; do
  case "$count" in ''|*[!0-9]*|0) echo "REFUSE: runtime thread counts must be positive integers" >&2; exit 2 ;; esac
  if [ "$count" -gt "$per_rank_threads" ]; then per_rank_threads=$count; fi
done
if [ $((np * per_rank_threads)) -gt "$runner_thread_budget" ]; then
  echo "REFUSE: MPI ranks ($np) x per-rank threads ($per_rank_threads) exceeds rank/thread budget $runner_thread_budget" >&2
  exit 2
fi

case "$bin" in
  plain)
    exe=$plain_exe; sparta_tag=${sparta_tag:-27Aug2026}
    sparta_commit=${sparta_commit:-95b9abaa8bd548991cc3c3f1c58b34722f7ade74} ;;
  *)
    exe=$bin; sparta_tag=${sparta_tag:-unknown}; sparta_commit=${sparta_commit:-unknown} ;;
esac

d=$run_root/$sub
if ! mkdir -p "$run_root" "$log_root"; then echo "REFUSE: cannot create run/log roots" >&2; exit 3; fi
say() { echo "$(date -u +%FT%TZ) $sub $*" | tee -a "$log_root/run_lean.log" >&2; }

if [ -e "$d" ]; then say "REFUSE: $d already exists"; exit 3; fi
if pgrep -x 'spa_.*' > /dev/null; then
  say "REFUSE: a SPARTA job is running: $(pgrep -ax 'spa_.*' | head -1 | cut -c1-120)"; exit 4
fi
if [ ! -x "$exe" ]; then say "REFUSE: solver is not executable: $exe"; exit 5; fi
if ! command -v "$python_cmd" >/dev/null 2>&1; then say "REFUSE: provenance writer is unavailable: $python_cmd"; exit 5; fi
if [ "$np" -gt 1 ] && ! command -v mpirun >/dev/null 2>&1; then say "REFUSE: mpirun is required when NPROC > 1"; exit 5; fi
repo_commit=$(git -C "$repo" rev-parse HEAD 2>/dev/null || true)
if ! printf '%s\n' "$repo_commit" | grep -Eq '^[0-9a-f]{40}$'; then
  say "REFUSE: Git provenance is unavailable or blank for $repo"; exit 6
fi

inputs=("$@" data/he.species data/he.vss)
sources=(); staged_names=(); roles=()
for f in "${inputs[@]}"; do
  src="$repo/cases/$casedir/$f"; [ -e "$src" ] || src="$repo/cases/$f"
  if [ ! -f "$src" ] || [ ! -r "$src" ]; then say "REFUSE: staged input is missing or unreadable: $f"; exit 7; fi
  name=$(basename "$f")
  for existing in "${staged_names[@]:-}"; do
    if [ "$existing" = "$name" ]; then say "REFUSE: staged inputs share the basename $name"; exit 7; fi
  done
  sources+=("$src"); staged_names+=("$name")
  case "$name" in
    "$deck") roles+=("deck") ;;
    *.species|*.vss) roles+=("collision-data") ;;
    *.surf) roles+=("geometry-or-station") ;;
    *) roles+=("case-input") ;;
  esac
done

varflags=()
IFS=',' read -ra kvs <<< "$vars"
for kv in "${kvs[@]}"; do
  if [[ "$kv" != *=* ]] || [ -z "${kv%%=*}" ]; then say "REFUSE: bad variable assignment: $kv"; exit 7; fi
  varflags+=(-var "${kv%%=*}" "${kv#*=}")
done
command_argv=("$exe" -in "$deck" "${varflags[@]}")

say "START np=$np deck=$deck exe=$exe"
if ! mkdir "$d"; then say "REFUSE: cannot reserve run directory: $d"; exit 3; fi
started=$(date -u +%FT%TZ)
if ! printf '%s\n' "$started" > "$d/start.stamp"; then say "REFUSE: cannot write start stamp"; exit 7; fi
for i in "${!sources[@]}"; do
  if ! tr -d '\r' < "${sources[$i]}" > "$d/${staged_names[$i]}"; then
    say "REFUSE: failed while staging ${sources[$i]}"; exit 7
  fi
done

if ! {
  echo "binary $exe"; sha256sum "$exe"
  echo "deck $deck"; echo "vars $vars"; echo "nproc $np"
  echo "inputs ${inputs[*]}"; echo "repo_commit $repo_commit"
} > "$d/manifest.txt"; then say "REFUSE: cannot write legacy receipt"; exit 8; fi

receipt=("$python_cmd" "$manifest_tool" start
  --run-dir "$d" --run-id "$sub" --case-id "$casedir" --repo-root "$repo"
  --solver-selector "$bin" --solver-tag "$sparta_tag" --solver-commit "$sparta_commit"
  --solver-exe "$exe" --started "$started" --mpi-ranks "$np"
  --runtime-control "OMP_NUM_THREADS=$OMP_NUM_THREADS"
  --runtime-control "OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS"
  --runtime-control "MKL_NUM_THREADS=$MKL_NUM_THREADS"
  --runtime-control "NUMEXPR_NUM_THREADS=$NUMEXPR_NUM_THREADS"
  --runtime-control "DSMC_THREAD_BUDGET=$runner_thread_budget")
if [ -z "$solver_receipt" ] && [ -f "$exe.build.json" ]; then solver_receipt=$exe.build.json; fi
if [ -n "$solver_receipt" ]; then receipt+=(--solver-receipt "$solver_receipt"); fi
if [ "$np" -eq 1 ]; then
  for arg in "${command_argv[@]}"; do receipt+=("--command-arg=$arg"); done
else
  receipt+=(--command-arg=mpirun)
  receipt+=(--command-arg=-np "--command-arg=$np")
  for arg in "${command_argv[@]}"; do receipt+=("--command-arg=$arg"); done
fi
for i in "${!sources[@]}"; do
  receipt+=(--input-record "${sources[$i]}"$'\t'"${staged_names[$i]}"$'\t'"${roles[$i]}")
done
if ! "${receipt[@]}"; then say "REFUSE: could not write provenance receipt"; exit 8; fi

cd "$d" || exit 8
if [ "$np" -eq 1 ]; then
  "${command_argv[@]}" &> run.log
else
  mpirun -np "$np" "${command_argv[@]}" &> run.log
fi
rc=$?
printf '%s\n' "$rc" > rc.sentinel
finished=$(date -u +%FT%TZ)
if ! "$python_cmd" "$manifest_tool" finish --run-dir "$d" --exit-code "$rc" --finished "$finished"; then
  say "ERROR: solver rc=$rc but final provenance update failed"; exit 8
fi
say "DONE rc=$rc"
echo "SENTINEL_RC=$rc RUN_DIR=$d"
echo "===TAIL"; tail -12 run.log
exit "$rc"
