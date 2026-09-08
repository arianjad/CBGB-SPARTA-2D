#!/usr/bin/env bash
# Install the tracer's pinned Julia release into this clone, without sudo.
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
root="$repo/.local"
[ "$(uname -s)" = Linux ] || { echo 'This installer supports Linux/WSL.' >&2; exit 2; }
case "$(uname -m)" in
  x86_64) arch=x86_64; urlarch=x64; sha=07d20c4c2518833e2265ca0acee15b355463361aa4efdab858dad826cf94325c ;;
  aarch64|arm64) arch=aarch64; urlarch=aarch64; sha=541d0c5a9378f8d2fc384bb8595fc6ffe20d61054629a6e314fb2f8dfe2f2ade ;;
  *) echo 'Unsupported architecture; use https://julialang.org/downloads/oldreleases/' >&2; exit 2 ;;
esac
[ ! -e "$root/julia-1.9.4" ] || { echo 'Julia directory already exists; use it or choose a fresh clone.' >&2; exit 2; }
mkdir -p "$root"
archive="$root/julia-1.9.4-linux-$arch.tar.gz"
curl --fail --location --retry 2 "https://julialang-s3.julialang.org/bin/linux/$urlarch/1.9/julia-1.9.4-linux-$arch.tar.gz" -o "$archive"
printf '%s  %s\n' "$sha" "$archive" | sha256sum --check -
tar -xzf "$archive" -C "$root"
"$root/julia-1.9.4/bin/julia" --version
