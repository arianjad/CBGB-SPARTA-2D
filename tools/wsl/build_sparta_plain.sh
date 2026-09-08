#!/usr/bin/env bash
# Build the pinned plain MPI SPARTA executable used by the 2D workflow.
set -euo pipefail

SPARTA_URL=${SPARTA_URL:-https://github.com/sparta/sparta.git}
SPARTA_TAG=${SPARTA_TAG:-27Aug2026}
SPARTA_COMMIT=${SPARTA_COMMIT:-95b9abaa8bd548991cc3c3f1c58b34722f7ade74}
SPARTA_SRC_ROOT=${SPARTA_SRC_ROOT:-$HOME/src/sparta-$SPARTA_TAG}
SPARTA_BUILD_ROOT=${SPARTA_BUILD_ROOT:-$HOME/build/sparta-$SPARTA_TAG-plain}
SPARTA_INSTALL_ROOT=${SPARTA_INSTALL_ROOT:-$HOME/opt/sparta-$SPARTA_TAG}
BUILD_JOBS=${BUILD_JOBS:-4}
TEST_JOBS=${TEST_JOBS:-1}

print_plan() {
  cat <<EOF
git clone $SPARTA_URL $SPARTA_SRC_ROOT
git -C $SPARTA_SRC_ROOT checkout --detach $SPARTA_COMMIT  # tag $SPARTA_TAG
cmake -S $SPARTA_SRC_ROOT/cmake -B $SPARTA_BUILD_ROOT -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$SPARTA_INSTALL_ROOT -DSPARTA_MACHINE=mpi -DSPARTA_ENABLE_TESTING=ON
cmake --build $SPARTA_BUILD_ROOT --parallel $BUILD_JOBS
ctest --test-dir $SPARTA_BUILD_ROOT --output-on-failure -j $TEST_JOBS
cmake --install $SPARTA_BUILD_ROOT
test -x $SPARTA_INSTALL_ROOT/bin/spa_mpi
EOF
}

if [ "${1:-}" = "--print-plan" ]; then print_plan; exit 0; fi
if [ "$#" -ne 0 ]; then echo "usage: $0 [--print-plan]" >&2; exit 2; fi

for command in git cmake ctest python3; do
  command -v "$command" >/dev/null || { echo "missing prerequisite: $command" >&2; exit 3; }
done
if [ -e "$SPARTA_BUILD_ROOT" ] || [ -e "$SPARTA_INSTALL_ROOT" ]; then
  echo "refusing existing build/install path; choose fresh roots" >&2
  exit 4
fi
if [ ! -e "$SPARTA_SRC_ROOT" ]; then
  git clone "$SPARTA_URL" "$SPARTA_SRC_ROOT"
elif [ ! -d "$SPARTA_SRC_ROOT/.git" ]; then
  echo "source root exists and is not a Git checkout: $SPARTA_SRC_ROOT" >&2
  exit 4
fi
actual=$(git -C "$SPARTA_SRC_ROOT" rev-parse HEAD 2>/dev/null || true)
if [ "$actual" != "$SPARTA_COMMIT" ]; then
  if [ -n "$(git -C "$SPARTA_SRC_ROOT" status --porcelain)" ]; then
    echo "refusing to change a dirty SPARTA source checkout: $SPARTA_SRC_ROOT" >&2
    exit 4
  fi
  git -C "$SPARTA_SRC_ROOT" fetch --tags origin
  git -C "$SPARTA_SRC_ROOT" checkout --detach "$SPARTA_COMMIT"
fi
actual=$(git -C "$SPARTA_SRC_ROOT" rev-parse HEAD)
[ "$actual" = "$SPARTA_COMMIT" ] || { echo "SPARTA commit mismatch: $actual" >&2; exit 5; }

cmake -S "$SPARTA_SRC_ROOT/cmake" -B "$SPARTA_BUILD_ROOT" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$SPARTA_INSTALL_ROOT" \
  -DSPARTA_MACHINE=mpi \
  -DSPARTA_ENABLE_TESTING=ON
cmake --build "$SPARTA_BUILD_ROOT" --parallel "$BUILD_JOBS"
ctest --test-dir "$SPARTA_BUILD_ROOT" --output-on-failure -j "$TEST_JOBS"
cmake --install "$SPARTA_BUILD_ROOT"
test -x "$SPARTA_INSTALL_ROOT/bin/spa_mpi"
solver_hash=$(sha256sum "$SPARTA_INSTALL_ROOT/bin/spa_mpi" | awk '{print $1}')
python3 - "$SPARTA_INSTALL_ROOT/bin/spa_mpi.build.json" "$SPARTA_URL" "$SPARTA_TAG" "$SPARTA_COMMIT" "$solver_hash" <<'PY'
import json, pathlib, sys
out, url, tag, commit, digest = sys.argv[1:]
pathlib.Path(out).write_text(json.dumps({
    "source_url": url, "tag": tag, "commit": commit,
    "executable_sha256": digest,
}, indent=2) + "\n", encoding="utf-8")
PY
printf '%s  %s\n' "$solver_hash" "$SPARTA_INSTALL_ROOT/bin/spa_mpi"
