#!/usr/bin/env bash
# Isolated RED/GREEN harness for run_lean.sh. It never starts SPARTA or WSL:
# each case invokes the tracked runner with scratch roots supplied through its
# documented environment, then uses an executable fake solver.
#
# Run:
#   bash cases/b5-lean/test_run_lean.sh
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
OWNED_TMP_ROOT=${TMPDIR:-/tmp}/dsmc-run-lean-harness
mkdir -p "$OWNED_TMP_ROOT"
TMP=$(mktemp -d "$OWNED_TMP_ROOT/case.XXXXXX")
TMP_CANON=$(cd "$TMP" && pwd -P)
OWNED_TMP_ROOT_CANON=$(cd "$OWNED_TMP_ROOT" && pwd -P)
case "$TMP_CANON" in
  "$OWNED_TMP_ROOT_CANON"/case.*) ;;
  *) echo "refusing unsafe temporary path: $TMP_CANON" >&2; exit 2 ;;
esac
cleanup() {
  case "$TMP_CANON" in
    "$OWNED_TMP_ROOT_CANON"/case.*) rm -rf -- "$TMP_CANON" ;;
    *) echo "refusing unsafe cleanup target: $TMP_CANON" >&2 ;;
  esac
}
trap cleanup EXIT
REPO="$TMP_CANON/repo"
RUNROOT="$TMP_CANON/runs"
LOGROOT="$TMP_CANON/logs"
FAKE="$TMP_CANON/fake-sparta"
FAKEBIN="$TMP_CANON/bin"
REAL_GIT=$(command -v git)

failures=0
check() {
  if ! "$@"; then
    echo "FAIL: $*" >&2
    failures=$((failures + 1))
  fi
}

cat > "$FAKE" <<'EOF'
#!/usr/bin/env bash
printf '%s\n' "$*" >> "${FAKE_LOG:?}"
if [ -n "${FAKE_BARRIER:-}" ]; then
  touch "$FAKE_BARRIER/$(basename "$PWD")"
  for attempt in $(seq 1 100); do
    [ "$(find "$FAKE_BARRIER" -type f | wc -l)" -ge 2 ] && break
    sleep 0.05
  done
  [ "$(find "$FAKE_BARRIER" -type f | wc -l)" -ge 2 ] || exit 99
fi
exit "${FAKE_RC:-0}"
EOF
chmod +x "$FAKE"

fresh() {
  for target in "$REPO" "$RUNROOT" "$LOGROOT"; do
    case "$target" in
      "$TMP_CANON"/repo|"$TMP_CANON"/runs|"$TMP_CANON"/logs) ;;
      *) echo "refusing unsafe fixture reset target: $target" >&2; exit 2 ;;
    esac
  done
  rm -rf -- "$REPO" "$RUNROOT" "$LOGROOT"
  mkdir -p "$REPO/cases/toy/data" "$RUNROOT" "$LOGROOT" "$FAKEBIN"
  printf 'deck\n' > "$REPO/cases/toy/deck.in"
  printf 'wall\r\n' > "$REPO/cases/toy/wall.surf"
  printf 'species\n' > "$REPO/cases/toy/data/he.species"
  printf 'vss\n' > "$REPO/cases/toy/data/he.vss"
  # Each case starts from an executable fake solver.  The noexe case may
  # remove that permission without contaminating the cases that follow.
  chmod +x "$FAKE"
  cat > "$FAKEBIN/git" <<EOF
#!/usr/bin/env bash
if [ "\${GIT_MODE:-real}" = unavailable ]; then exit 127; fi
if [ "\${GIT_MODE:-real}" = blank ]; then exit 0; fi
exec "$REAL_GIT" "\$@"
EOF
  chmod +x "$FAKEBIN/git"
  # Another visible SPARTA process must not veto an independent run.
  printf '#!/usr/bin/env bash\necho "123 spa_mpi"\nexit 0\n' > "$FAKEBIN/pgrep"
  chmod +x "$FAKEBIN/pgrep"
  : > "$TMP/fake.log"
}

make_git_repo() {
  git -C "$REPO" init -q
  git -C "$REPO" config user.email runner-test@example.invalid
  git -C "$REPO" config user.name runner-test
  git -C "$REPO" add cases
  git -C "$REPO" commit -qm fixture
}

run() {
  local sub=$1; shift
  GIT_MODE="${GIT_MODE:-real}" PATH="$FAKEBIN:$PATH" \
    FAKE_LOG="$TMP/fake.log" FAKE_RC="${FAKE_RC:-0}" \
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
    DSMC_REPO_ROOT="$REPO" DSMC_RUN_ROOT="$RUNROOT" DSMC_LOG_ROOT="$LOGROOT" \
    DSMC_MANIFEST_TOOL="$ROOT/tools/lean_manifest.py" \
    SPARTA_TAG="test-tag" SPARTA_COMMIT="1111111111111111111111111111111111111111" \
    bash "$ROOT/cases/b5-lean/run_lean.sh" "$sub" "${RUN_NP:-1}" deck.in "A=1" "$FAKE" toy "$@"
}

expect_refusal() {
  local name=$1; shift
  if "$@"; then
    echo "FAIL: $name launched" >&2
    failures=$((failures + 1))
  fi
  check test ! -s "$TMP/fake.log"
}

# Existing green behavior: a committed repository produces a nonblank receipt,
# stages every named input, and returns the solver's zero status.
fresh; make_git_repo; FAKE_RC=0 run valid deck.in wall.surf
check test "$(cat "$RUNROOT/valid/rc.sentinel")" = 0
check grep -Eq '^repo_commit [0-9a-f]{40}$' "$RUNROOT/valid/manifest.txt"
check test -s "$RUNROOT/valid/deck.in"
check test -s "$RUNROOT/valid/he.species"
check test -s "$RUNROOT/valid/manifest.json"
python3 - "$RUNROOT/valid/manifest.json" "$REPO" "$FAKE" <<'PY'
import hashlib, json, pathlib, sys
m = json.loads(pathlib.Path(sys.argv[1]).read_text())
root, solver = pathlib.Path(sys.argv[2]), pathlib.Path(sys.argv[3])
assert m["status"] == "complete" and m["failure"] is None
assert m["case_id"] == "toy" and m["repository"]["dirty"] is False
assert m["solver"]["executable"] == str(solver)
assert m["solver"]["tag"] == "test-tag"
assert m["solver"]["commit"] == "1" * 40
assert m["solver"]["source_identity"] == "declared"
assert m["runtime_controls"]["mpi_ranks"] == 1
assert m["runtime_controls"]["OMP_NUM_THREADS"] == "1"
assert m["runtime_controls"]["DSMC_THREAD_BUDGET"] == "16"
assert m["command"][-5:] == ["-in", "deck.in", "-var", "A", "1"]
by_name = {x["staged_path"]: x for x in m["inputs"]}
for name in ("deck.in", "wall.surf", "he.species", "he.vss"):
    item = by_name[name]
    source = root / item["source_path"]
    staged = pathlib.Path(sys.argv[1]).parent / name
    assert item["source_sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert item["staged_sha256"] == hashlib.sha256(staged.read_bytes()).hexdigest()
assert by_name["wall.surf"]["source_sha256"] != by_name["wall.surf"]["staged_sha256"]
assert any(x["path"] == "run.log" for x in m["outputs"])
assert not set(by_name).intersection(x["path"] for x in m["outputs"])
PY

# Both independent solvers must be active before either can finish.
fresh; make_git_repo; mkdir "$TMP/barrier"
export FAKE_BARRIER="$TMP/barrier"
run sweep-a deck.in wall.surf & first=$!
run sweep-b deck.in wall.surf & second=$!
check wait "$first"
check wait "$second"
unset FAKE_BARRIER
check test "$(cat "$RUNROOT/sweep-a/rc.sentinel")" = 0
check test "$(cat "$RUNROOT/sweep-b/rc.sentinel")" = 0

# Dirty source state is recorded rather than silently attributed to HEAD.
fresh; make_git_repo; printf 'untracked\n' > "$REPO/local-note.txt"
FAKE_RC=0 run dirty deck.in wall.surf
python3 - "$RUNROOT/dirty/manifest.json" <<'PY'
import json, pathlib, sys
assert json.loads(pathlib.Path(sys.argv[1]).read_text())["repository"]["dirty"] is True
PY

# Attempt every refusal case even after an earlier guard fails.
fresh; make_git_repo; mkdir "$RUNROOT/partial"; touch "$RUNROOT/partial/old-output"
expect_refusal "partial output directory" run partial deck.in wall.surf
check test -f "$RUNROOT/partial/old-output"

# A completed directory is equally immutable evidence.
fresh; make_git_repo; mkdir "$RUNROOT/completed"; printf '0\n' > "$RUNROOT/completed/rc.sentinel"
expect_refusal "completed output directory" run completed deck.in wall.surf
check test "$(cat "$RUNROOT/completed/rc.sentinel")" = 0

# A missing staged file must refuse before the fake solver is called.
fresh; make_git_repo
expect_refusal "missing staged input" run missing deck.in absent.surf
check test ! -e "$RUNROOT/missing"

# A non-executable solver must refuse before staging a launch receipt.
fresh; make_git_repo; chmod -x "$FAKE"
if [ -x "$FAKE" ]; then
  echo "SKIP: Git Bash does not expose chmod -x for this fixture"
else
  expect_refusal "non-executable solver" run noexe deck.in wall.surf
  check test ! -e "$RUNROOT/noexe"
fi

# Both an unavailable Git repository and a successful Git command with blank
# output must refuse; the latter catches error suppression that masks a bad
# provenance receipt.
fresh; make_git_repo
GIT_MODE=unavailable expect_refusal "unavailable Git provenance" run nogit deck.in wall.surf
check test ! -e "$RUNROOT/nogit"

fresh; make_git_repo
GIT_MODE=blank expect_refusal "blank Git provenance" run blankgit deck.in wall.surf
check test ! -e "$RUNROOT/blankgit"

# The configurable rank/thread budget defaults to the managed DSMC cap of 16.
fresh; make_git_repo
if RUN_NP=17 run overbudget deck.in wall.surf > "$TMP/overbudget.out" 2>&1; then
  echo "FAIL: over-budget launch returned success" >&2; failures=$((failures + 1))
fi
check grep -q "rank/thread budget" "$TMP/overbudget.out"
check test ! -e "$RUNROOT/overbudget"

# A fake solver failure must survive the runner's status/log tailing path.
fresh; make_git_repo
rc=0
if FAKE_RC=7 run solverfail deck.in wall.surf; then
  echo "FAIL: solver failure returned success" >&2
  failures=$((failures + 1))
else
  rc=$?
fi
check test "$rc" = 7
check test "$(cat "$RUNROOT/solverfail/rc.sentinel")" = 7
python3 - "$RUNROOT/solverfail/manifest.json" <<'PY'
import json, pathlib, sys
m = json.loads(pathlib.Path(sys.argv[1]).read_text())
assert m["status"] == "failed"
assert m["failure"]["exit_code"] == 7
assert m["finished_utc"]
PY

if [ "$failures" -ne 0 ]; then
  echo "FAIL run_lean isolated harness: $failures check(s) failed" >&2
  exit 1
fi
echo "PASS run_lean isolated harness"
