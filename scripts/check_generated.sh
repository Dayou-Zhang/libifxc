#!/bin/sh
# SPDX-License-Identifier: MPL-2.0
set -eu

script_dir=$(CDPATH= cd -- "$(dirname "$0")" && pwd)
repo_root=$(CDPATH= cd -- "${script_dir}/.." && pwd)
python_cmd=${1:-python3}

"$python_cmd" "$repo_root/scripts/check_ml25_manifest.py" \
  --source "$repo_root/maple/if_mgga/mgga_xc_ml25.mpl" \
  --feature-def "$repo_root/src/features/ifxc_ml25_features.def"

"$python_cmd" "$repo_root/scripts/check_ml26_manifest.py" \
  --source "$repo_root/maple/if_mgga/mgga_xc_ml26.mpl" \
  --feature-def "$repo_root/src/features/ifxc_ml26_features.def"
