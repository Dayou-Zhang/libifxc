#!/bin/sh
set -eu

script_dir=$(CDPATH= cd -- "$(dirname "$0")" && pwd)
repo_root=$(CDPATH= cd -- "${script_dir}/.." && pwd)
build_dir="${1:-${repo_root}/build}"

cmake --build "$build_dir" --target generate

cd "$repo_root"
git diff --exit-code -- maple src/features src/generated
