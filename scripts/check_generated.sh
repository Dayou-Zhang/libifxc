#!/bin/sh
set -eu

script_dir=$(CDPATH= cd -- "$(dirname "$0")" && pwd)
repo_root=$(CDPATH= cd -- "${script_dir}/.." && pwd)
build_dir="${1:-${repo_root}/build}"
tmpdir=$(mktemp -d)

cleanup() {
  rm -rf "$tmpdir"
}
trap cleanup EXIT INT TERM

cp -R "$repo_root/src/generated" "$tmpdir/generated.before"
cmake --build "$build_dir" --target generate

diff -ru "$tmpdir/generated.before" "$repo_root/src/generated"
