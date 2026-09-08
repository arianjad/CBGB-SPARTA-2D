#!/usr/bin/env bash
# Convert one completed 2D helium run, trace the configured heavy species,
# and emit exit/crossing metrics plus helium/molecule figures.
#
# Usage: tools/run_2d_standalone.sh HELIUM_RUN OUT_DIR [SOURCE_CONFIG]
# Runtime controls: N (1000), THREADS (4), SEED (42), PYTHON (python3),
# JULIA (julia), TRACER_THREAD_BUDGET (4), HEATBINS (8,4), TRAJPRINT (0).
set -euo pipefail

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
  echo "usage: $0 HELIUM_RUN OUT_DIR [SOURCE_CONFIG]" >&2
  exit 2
fi

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
helium_run=$1
out=$2
config=${3:-$repo/cases/b5-lean/baf-hot-source.conf}
python_cmd=${PYTHON:-python3}
julia_cmd=${JULIA:-julia}
n=${N:-1000}
threads=${THREADS:-4}
seed=${SEED:-42}
thread_budget=${TRACER_THREAD_BUDGET:-4}
heatbins=${HEATBINS:-8,4}
trajprint=${TRAJPRINT:-0}

case "$n" in ''|*[!0-9]*|0) echo "REFUSE: N must be a positive integer" >&2; exit 2 ;; esac
case "$threads" in ''|*[!0-9]*|0) echo "REFUSE: THREADS must be a positive integer" >&2; exit 2 ;; esac
case "$thread_budget" in ''|*[!0-9]*|0) echo "REFUSE: TRACER_THREAD_BUDGET must be a positive integer" >&2; exit 2 ;; esac
case "$seed" in ''|*[!0-9]*) echo "REFUSE: SEED must be a non-negative integer" >&2; exit 2 ;; esac
case "$trajprint" in ''|*[!0-9]*) echo "REFUSE: TRAJPRINT must be a non-negative integer" >&2; exit 2 ;; esac
if [ "$threads" -gt "$thread_budget" ]; then
  echo "REFUSE: THREADS=$threads exceeds TRACER_THREAD_BUDGET=$thread_budget" >&2
  exit 2
fi
case "$heatbins" in *','*) ;; *) echo "REFUSE: HEATBINS must be NZ,NR" >&2; exit 2 ;; esac

[ -d "$helium_run" ] || { echo "missing helium run directory: $helium_run" >&2; exit 3; }
[ -f "$helium_run/field.grid" ] || { echo "missing helium field: $helium_run/field.grid" >&2; exit 3; }
compgen -G "$helium_run/cell*.surf" >/dev/null || {
  echo "missing staged cell*.surf in helium run: $helium_run" >&2; exit 3;
}
[ -f "$config" ] || { echo "missing source config: $config" >&2; exit 3; }
[ ! -e "$out" ] || { echo "REFUSE: output path already exists: $out" >&2; exit 4; }
command -v "$python_cmd" >/dev/null || { echo "missing Python command: $python_cmd" >&2; exit 3; }
command -v "$julia_cmd" >/dev/null || { echo "missing Julia command: $julia_cmd" >&2; exit 3; }

# A populated field can exist while its solver is still running or has failed.
"$python_cmd" - "$helium_run" <<'PY'
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
try:
    status = json.loads((root / 'manifest.json').read_text())['status']
    rc = (root / 'rc.sentinel').read_text().strip()
except (OSError, ValueError, KeyError) as exc:
    raise SystemExit(f'incomplete helium receipt: {exc}')
if status != 'complete' or rc != '0':
    raise SystemExit(f'helium run is not complete: status={status}, rc={rc}')
PY

helium_run=$(cd "$helium_run" && pwd -P)
config=$(cd "$(dirname "$config")" && pwd -P)/$(basename "$config")
config_origin=$config
config_text=$(cat "$config_origin")

# Parse one in-memory snapshot and write that same snapshot into the receipt,
# so a concurrent edit of the repository copy cannot change the values used.
# This configuration is data, not a shell script.
while IFS='=' read -r key value; do
  key=${key//$'\r'/}; value=${value//$'\r'/}
  case "$key" in ''|'#'*) continue ;; esac
  case "$key" in
    MASS_U|CROSS_SECTION_M2|SOURCE_T_K|SPAWN_MODE|SPAWN_SIZE_M|SPAWN_R_M|SPAWN_Z_M|SPAWN_CLIP|SAMPLER|OBS_X_M)
      printf -v "$key" '%s' "$value" ;;
    *) echo "REFUSE: unknown source-config key: $key" >&2; exit 2 ;;
  esac
done <<< "$config_text"
for key in MASS_U CROSS_SECTION_M2 SOURCE_T_K SPAWN_MODE SPAWN_SIZE_M SPAWN_R_M SPAWN_Z_M SPAWN_CLIP SAMPLER OBS_X_M; do
  [ -n "${!key:-}" ] || { echo "REFUSE: source config is missing $key" >&2; exit 2; }
done
case "$SPAWN_MODE" in point|gaussball|uniformball) ;; *)
  echo "REFUSE: unsupported SPAWN_MODE=$SPAWN_MODE" >&2; exit 2 ;;
esac
case "$SAMPLER" in exact|table) ;; *) echo "REFUSE: unsupported SAMPLER=$SAMPLER" >&2; exit 2 ;; esac
case "$SPAWN_CLIP" in 0|1) ;; *) echo "REFUSE: SPAWN_CLIP must be 0 or 1" >&2; exit 2 ;; esac

mkdir -p "$out"
out=$(cd "$out" && pwd -P)
mkdir -p "$out/field" "$out/trace" "$out/figures"
printf '%s\n' "$config_text" > "$out/source.conf"
config="$out/source.conf"

repo_commit=$(git -C "$repo" rev-parse HEAD 2>/dev/null || true)
repo_dirty=false
[ -z "$(git -C "$repo" status --porcelain --untracked-files=all 2>/dev/null)" ] || repo_dirty=true
julia_version=$("$julia_cmd" --version 2>&1 | head -1)
python_version=$("$python_cmd" --version 2>&1 | head -1)
created=$(date -u +%FT%TZ)
manifest="$out/manifest.json"
stage=initializing

export STANDALONE_MANIFEST="$manifest" STANDALONE_REPO="$repo" STANDALONE_REPO_COMMIT="$repo_commit"
export STANDALONE_REPO_DIRTY="$repo_dirty" STANDALONE_HELIUM_RUN="$helium_run" STANDALONE_CONFIG="$config"
export STANDALONE_CONFIG_ORIGIN="$config_origin"
export STANDALONE_N="$n" STANDALONE_THREADS="$threads" STANDALONE_THREAD_BUDGET="$thread_budget"
export STANDALONE_SEED="$seed" STANDALONE_TRAJPRINT="$trajprint" STANDALONE_OBS_X="$OBS_X_M" STANDALONE_JULIA="$julia_cmd"
export STANDALONE_JULIA_VERSION="$julia_version" STANDALONE_PYTHON="$python_cmd"
export STANDALONE_PYTHON_VERSION="$python_version" STANDALONE_CREATED="$created" STANDALONE_OUT="$out"
export MASS_U CROSS_SECTION_M2 SOURCE_T_K SPAWN_MODE SPAWN_SIZE_M SPAWN_R_M SPAWN_Z_M SPAWN_CLIP SAMPLER

write_manifest() {
  local status=$1 current=$2 exit_code=${3:-}
  STANDALONE_STATUS="$status" STANDALONE_STAGE="$current" STANDALONE_EXIT="$exit_code" \
    "$python_cmd" - "$manifest" "$repo" "$helium_run" "$config" "$out" <<'PY'
import hashlib, json, os
import sys
from datetime import datetime, timezone
from pathlib import Path

p, repo_arg, helium_arg, config_arg, out_arg = map(Path, sys.argv[1:6])

def digest(q):
    q = Path(q)
    if not q.is_file():
        return None
    h = hashlib.sha256()
    with q.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

if p.exists():
    data = json.loads(p.read_text(encoding="utf-8"))
else:
    repo = repo_arg
    he = helium_arg
    cfg = config_arg
    inputs = [he / "field.grid", he / "manifest.json", cfg,
              repo / "tools/run_2d_standalone.sh",
              repo / "tools/field2tracer.py", repo / "tools/plot_fields.py",
              repo / "tools/plot_fields_b5.py", repo / "tools/plot_trajectories.py",
              repo / "tracer/accumulators/crossing.jl",
              repo / "tracer/ParticleTracing.jl", repo / "tracer/Project.toml",
              repo / "tracer/Manifest.toml", repo / "tools/tracer_analyze.py",
              repo / "tools/tracer_common.py"]
    data = {
        "schema": "dsmc-2d-standalone-v1",
        "created_utc": os.environ["STANDALONE_CREATED"],
        "repo_root": str(repo),
        "repo_commit": os.environ["STANDALONE_REPO_COMMIT"] or None,
        "repo_dirty": os.environ["STANDALONE_REPO_DIRTY"] == "true",
        "helium_run": str(he),
        "helium_manifest_status": None,
        "source_config": str(cfg),
        "source_config_origin": os.environ["STANDALONE_CONFIG_ORIGIN"],
        "n_particles": int(os.environ["STANDALONE_N"]),
        "particle_seed": int(os.environ["STANDALONE_SEED"]),
        "trajectory_print_particles": int(os.environ["STANDALONE_TRAJPRINT"]),
        "julia_threads": int(os.environ["STANDALONE_THREADS"]),
        "thread_budget": int(os.environ["STANDALONE_THREAD_BUDGET"]),
        "observation_plane_m": float(os.environ["STANDALONE_OBS_X"]),
        "molecular_source": {
            "mass_u": float(os.environ["MASS_U"]),
            "cross_section_m2": float(os.environ["CROSS_SECTION_M2"]),
            "temperature_K": float(os.environ["SOURCE_T_K"]),
            "spawn_mode": os.environ["SPAWN_MODE"],
            "spawn_size_m": float(os.environ["SPAWN_SIZE_M"]),
            "spawn_r_m": float(os.environ["SPAWN_R_M"]),
            "spawn_z_m": float(os.environ["SPAWN_Z_M"]),
            "spawn_clip": int(os.environ["SPAWN_CLIP"]),
            "sampler": os.environ["SAMPLER"],
        },
        "julia": {"command": os.environ["STANDALONE_JULIA"],
                  "version": os.environ["STANDALONE_JULIA_VERSION"]},
        "python": {"command": os.environ["STANDALONE_PYTHON"],
                   "version": os.environ["STANDALONE_PYTHON_VERSION"]},
        "input_sha256": {str(q): digest(q) for q in inputs},
    }
    hm = he / "manifest.json"
    if hm.exists():
        try:
            data["helium_manifest_status"] = json.loads(hm.read_text(encoding="utf-8")).get("status")
        except Exception:
            data["helium_manifest_status"] = "unreadable"
data["status"] = os.environ["STANDALONE_STATUS"]
data["current_stage"] = os.environ["STANDALONE_STAGE"]
data["updated_utc"] = datetime.now(timezone.utc).isoformat()
if os.environ["STANDALONE_EXIT"]:
    data["exit_code"] = int(os.environ["STANDALONE_EXIT"])
    root = out_arg
    data["artifacts"] = sorted(str(q.relative_to(root)) for q in root.rglob("*")
                               if q.is_file() and q != p)
p.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
PY
}

finish() {
  rc=$?
  trap - EXIT
  printf '%s\n' "$rc" > "$out/rc.sentinel"
  if [ "$rc" -eq 0 ]; then
    write_manifest complete complete 0
  else
    write_manifest failed "$stage" "$rc" || true
  fi
  exit "$rc"
}
trap finish EXIT
write_manifest running "$stage"

record_command() {
  printf '%s: ' "$1" >> "$out/commands.sh"
  shift
  printf '%q ' "$@" >> "$out/commands.sh"
  printf '\n' >> "$out/commands.sh"
}

stage=convert_field
write_manifest running "$stage"
convert=("$python_cmd" "$repo/tools/field2tracer.py" "$helium_run" --out "$out/field")
record_command "$stage" "${convert[@]}"
"${convert[@]}" > "$out/convert.stdout" 2> "$out/convert.stderr"

stage=verify_field
write_manifest running "$stage"
"$python_cmd" - "$out/field/provenance.json" "$out/field/DS2FF.DAT" <<'PY'
import json, sys
from pathlib import Path

prov = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
step = prov.get("source_timestep")
if step is None or step <= 0:
    raise SystemExit("refusing tracer launch: final field snapshot is step 0; run at least one complete field-averaging interval")
populated = 0
with Path(sys.argv[2]).open(encoding="utf-8") as f:
    for i, line in enumerate(f):
        if i < 9:
            continue
        cols = line.split()
        if len(cols) >= 4 and float(cols[2]) > 0 and float(cols[3]) > 0:
            populated += 1
if populated < 100:
    raise SystemExit(f"refusing tracer launch: only {populated} populated field cells; tracer requires at least 100")
print(f"verified tracer field: timestep {step}, {populated} populated cells")
PY

stage=plot_helium
write_manifest running "$stage"
plot_he=("$python_cmd" "$repo/tools/plot_fields_b5.py" "$helium_run" --frac 0 --outdir "$out/figures")
record_command "$stage" "${plot_he[@]}"
"${plot_he[@]}" > "$out/helium-plot.stdout" 2> "$out/helium-plot.stderr"

stage=trace_crossings
write_manifest running "$stage"
trace=("$julia_cmd" --project="$repo/tracer" --threads="$threads"
       "$repo/tracer/accumulators/crossing.jl" "$out/trace" "$OBS_X_M"
       "$out/field/cell.surfs" "$out/field/DS2FF.DAT"
       -n "$n" --seed "$seed" -M "$MASS_U" --sigma "$CROSS_SECTION_M2" -T "$SOURCE_T_K"
       --spawn "$SPAWN_MODE" --spawnsize "$SPAWN_SIZE_M" -r "$SPAWN_R_M" -z "$SPAWN_Z_M"
       --spawnclip "$SPAWN_CLIP" --sampler "$SAMPLER" --saveall 1 --trajprint "$trajprint"
       --spawnout "$out/trace/spawn.csv")
record_command "$stage" "${trace[@]}"
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  "${trace[@]}" > "$out/trace/tracer.stdout" 2> "$out/trace/tracer.stderr"

if [ "$trajprint" -gt 0 ]; then
  stage=plot_trajectories
  write_manifest running "$stage"
  paths="$out/trace/paths.out"
  {
    awk -v n="$trajprint" '$1 < 0 && -$1 <= n {print}' "$out/trace/tracer.stdout"
    awk -v n="$trajprint" 'NR > 1 && $1 >= 1 && $1 <= n {print}' "$out/trace/legs.out"
  } > "$paths"
  plot_paths=("$python_cmd" "$repo/tools/plot_trajectories.py" "$paths" "$out/field/cell.surfs"
              --outdir "$out/figures" --label standalone-2d)
  record_command "$stage" "${plot_paths[@]}"
  "${plot_paths[@]}" > "$out/trajectory-plot.stdout" 2> "$out/trajectory-plot.stderr"
fi

stage=analyze_exit
write_manifest running "$stage"
analyze=("$python_cmd" "$repo/tools/tracer_analyze.py" "$out/trace/legs.out"
         "$out/field/cell.surfs" -N "$n" --outdir "$out/figures"
         --json "$out/exit.json" --label standalone-2d)
if [ -f "$out/trace/spawn.csv" ] && [ "$(($(wc -l < "$out/trace/spawn.csv") - 1))" -eq "$n" ]; then
  analyze+=(--spawnout "$out/trace/spawn.csv" --heatbins "$heatbins")
fi
record_command "$stage" "${analyze[@]}"
"${analyze[@]}" > "$out/exit.stdout" 2> "$out/exit.stderr"

stage=score_common_plane
write_manifest running "$stage"
common=("$python_cmd" "$repo/tools/tracer_common.py" "$out/trace" "$out/field/cell.surfs"
        --xobs "$OBS_X_M" --json "$out/common.json")
record_command "$stage" "${common[@]}"
"${common[@]}" > "$out/common.stdout" 2> "$out/common.stderr"

stage=complete
