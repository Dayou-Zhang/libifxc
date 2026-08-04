#!/usr/bin/env python3
# SPDX-License-Identifier: MPL-2.0
"""Validate ML25 manifest metadata without regenerating formula code."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


FEATURE_DEF_RE = re.compile(
    r'^IFXC_ML25_FEATURE\(\s*(\d+)\s*,\s*([A-Z0-9_]+)\s*,\s*"([^"]+)"\s*,\s*"([^"]+)"\s*,\s*(IFXC_FEATURE_KIND_[A-Z_]+)\s*\)\s*$'
)
INCLUDE_RE = re.compile(r'^\$include\s+"([^"]+)"\s*$')
KEY_STRING_RE = re.compile(r'^"([a-z][a-z0-9_]*)"[,]?\s*$')

REQUIRED_SUPPORT_FILES = [
    Path("maple/util.mpl"),
    Path("maple/lda_exc/lda_x.mpl"),
    Path("maple/lda_exc/lda_c_pw.mpl"),
    Path("maple/gga_exc/gga_c_lyp.mpl"),
    Path("maple/gga_exc/gga_c_pbe.mpl"),
    Path("maple/mgga_exc/mgga_x_mn12.mpl"),
    Path("maple/mgga_exc/mgga_c_m08.mpl"),
]
IF_MGGA_HELPERS = [
    Path("maple/if_mgga/ml25_lak.mpl"),
    Path("maple/if_mgga/ml25_lyp.mpl"),
    Path("maple/if_mgga/ml25_mn12_terms.mpl"),
    Path("maple/if_mgga/ml25_m08_terms.mpl"),
]

FORBIDDEN_GENERATED_METADATA = [
    Path("src/generated/ifxc_mgga_xc_ml25_features.c"),
    Path("src/generated/ifxc_mgga_xc_ml25_features.h"),
]
FORBIDDEN_ML25_MAPLE_PATTERNS = [
    "maple/mgga_xc_ml25.mpl",
    "maple/mgga_exc/mgga_x_lak.mpl",
    "maple/mgga_exc/mgga_c_lak.mpl",
    "maple/mgga_exc/mgga_xc_ml25_*.mpl",
]
FORBIDDEN_ML25_GENERATED_PATTERNS = [
    "src/ifxc_ml25_kernels.c",
    "src/maple2c/mgga_exc/mgga_x_lak.c",
    "src/maple2c/mgga_exc/mgga_c_lak.c",
    "src/maple2c/mgga_exc/mgga_xc_ml25*.c",
]
IF_MGGA_GENERATED_C = Path("src/maple2c/if_mgga/mgga_xc_ml25.c")


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def parse_feature_def(path: Path) -> list[tuple[int, str, str, str, str]]:
    features: list[tuple[int, str, str, str, str]] = []
    for line in read_text(path).splitlines():
        line = line.strip()
        if not line or line.startswith("/*"):
            continue
        match = FEATURE_DEF_RE.match(line)
        if match:
            features.append(
                (
                    int(match.group(1)),
                    match.group(2),
                    match.group(3),
                    match.group(4),
                    match.group(5),
                )
            )
    return features


def parse_maple_key_list(path: Path) -> list[tuple[int, str]]:
    keys: list[tuple[int, str]] = []
    in_vector = False
    for raw_line in read_text(path).splitlines():
        line = raw_line.strip()
        if line.startswith("ml25_feature_keys := ["):
            in_vector = True
            continue
        if not in_vector:
            continue
        if line.startswith("]:"):
            break
        if not line or line.startswith("#"):
            continue
        match = KEY_STRING_RE.match(line)
        if match:
            keys.append((len(keys), match.group(1)))
    return keys

def fail(message: str) -> int:
    print(message, file=sys.stderr)
    return 1


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True)
    parser.add_argument("--feature-def", required=True)
    args = parser.parse_args(argv)

    source_path = Path(args.source)
    feature_def_path = Path(args.feature_def)
    repo_root = feature_def_path.resolve().parents[2]

    source_text = read_text(source_path)
    if "ifxc_feature_set_key := \"ml25\":" not in source_text:
        return fail("Maple source does not advertise feature set ml25")
    if "ifxc_feature_count := 66:" not in source_text:
        return fail("Maple source does not advertise a 66-feature vector")
    if "(* type: if_mgga *)" not in source_text:
        return fail("ML25 source must be an IF-MGGA vector source")
    for directive in (
        "(* feature_set: ml25 *)",
        "(* nfeatures: 66 *)",
        "(* max_order: 3 *)",
        "(* variables: rho sigma tau *)",
    ):
        if directive not in source_text:
            return fail(f"missing IF-MGGA source directive: {directive}")

    feature_defs = parse_feature_def(feature_def_path)
    if len(feature_defs) != 66:
        return fail(f"expected 66 feature definitions, found {len(feature_defs)}")
    if [item[0] for item in feature_defs] != list(range(66)):
        return fail("feature definition indexes are not contiguous")

    maple_keys = parse_maple_key_list(source_path)
    if len(maple_keys) != 66:
        return fail(f"expected 66 Maple feature entries, found {len(maple_keys)}")
    if [item[0] for item in maple_keys] != list(range(66)):
        return fail("Maple feature indexes are not contiguous")

    feature_keys = [item[2] for item in feature_defs]
    maple_feature_keys = [item[1] for item in maple_keys]
    if feature_keys != maple_feature_keys:
        return fail("feature def keys and Maple feature keys differ")

    for support_path in REQUIRED_SUPPORT_FILES:
        if not (repo_root / support_path).exists():
            return fail(f"missing Maple support file: {support_path}")
    for helper_path in IF_MGGA_HELPERS:
        full_helper_path = repo_root / helper_path
        if not full_helper_path.exists():
            return fail(f"missing IF-MGGA helper file: {helper_path}")
        helper_include = f'$include "{helper_path.name}"'
        if helper_include not in source_text:
            return fail(f"IF-MGGA source does not include helper: {helper_path.name}")
        helper_text = read_text(full_helper_path)
        if re.search(r'(^|\n)\s*f\s*:=', helper_text):
            return fail(f"IF-MGGA helper assigns generic f: {helper_path}")
    for pattern in FORBIDDEN_ML25_MAPLE_PATTERNS:
        matches = sorted(repo_root.glob(pattern))
        if matches:
            return fail(f"legacy ML25 Maple wrapper remains: {matches[0].relative_to(repo_root)}")
    for pattern in FORBIDDEN_ML25_GENERATED_PATTERNS:
        matches = sorted(repo_root.glob(pattern))
        if matches:
            return fail(f"legacy ML25 generated file remains: {matches[0].relative_to(repo_root)}")
    for generated_path in FORBIDDEN_GENERATED_METADATA:
        if (repo_root / generated_path).exists():
            return fail(f"redundant generated metadata file remains: {generated_path}")

    combined_generated_c = IF_MGGA_GENERATED_C
    combined_path = repo_root / combined_generated_c
    if not combined_path.exists():
        return fail(f"missing combined generated ML25 file: {combined_generated_c}")
    combined_text = read_text(combined_path)
    for marker in (
        "This file was generated automatically with scripts/maple2c.py.",
        "Maple source      : maple/if_mgga/mgga_xc_ml25.mpl",
    ):
        if marker not in combined_text:
            return fail(f"{combined_generated_c} has a non-reproducible generated path: {marker}")
    for absolute_home in ("/home/", "/Users/", "\\\\Users\\\\"):
        if absolute_home in combined_text:
            return fail(f"{combined_generated_c} exposes an absolute home path")
    if "Error," in combined_text:
        return fail(f"Maple error text remains in {combined_generated_c}")
    if "Type of functional: if_mgga" not in combined_text:
        return fail(f"{combined_generated_c} was not generated as if_mgga")
    if "#define ifxc_maple2c_order 3" not in combined_text:
        return fail(f"{combined_generated_c} was not generated through third derivatives")
    if "#define maple2c_order" in combined_text:
        return fail(f"{combined_generated_c} defines unprefixed maple2c_order")
    for flag in ("IFXC_MGGA_FLAGS_HAVE_VXC", "IFXC_MGGA_FLAGS_HAVE_FXC", "IFXC_MGGA_FLAGS_HAVE_KXC"):
        if flag not in combined_text:
            return fail(f"{combined_generated_c} does not advertise {flag}")
    if "out->zk[ip*p->dim.zk + 65]" not in combined_text:
        return fail(f"{combined_generated_c} does not emit all 66 order-0 features")
    if "out->vlapl" in combined_text:
        return fail(f"{combined_generated_c} emits unsupported Laplacian derivatives")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
