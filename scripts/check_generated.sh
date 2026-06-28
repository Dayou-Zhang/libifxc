#!/bin/sh
set -eu

script_dir=$(CDPATH= cd -- "$(dirname "$0")" && pwd)
repo_root=$(CDPATH= cd -- "${script_dir}/.." && pwd)

python3 "$repo_root/scripts/check_ml25_manifest.py" \
  --source "$repo_root/maple/if_mgga/mgga_xc_ml25.mpl" \
  --feature-def "$repo_root/src/features/ifxc_ml25_features.def"
