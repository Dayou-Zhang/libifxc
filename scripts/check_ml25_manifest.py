#!/usr/bin/env python3
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
ASSIGN_RE = re.compile(r'^([a-z0-9_]+)\s*:=\s*f:\s*$')
MAPLE_KEY_RE = re.compile(r'#\s*(\d+)\s*,\s*(ml25\.[a-z0-9_]+)\s*$')

REQUIRED_SUPPORT_FILES = [
    Path("maple/util.mpl"),
    Path("maple/lda_exc/lda_x.mpl"),
    Path("maple/lda_exc/lda_c_pw.mpl"),
    Path("maple/gga_exc/gga_c_lyp.mpl"),
    Path("maple/gga_exc/gga_c_pbe.mpl"),
    Path("maple/mgga_exc/mgga_x_mn12.mpl"),
    Path("maple/mgga_exc/mgga_c_m08.mpl"),
]

FORBIDDEN_GENERATED_METADATA = [
    Path("src/generated/ifxc_mgga_xc_ml25_features.c"),
    Path("src/generated/ifxc_mgga_xc_ml25_features.h"),
]


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


def maple_include_for_key(feature_key: str) -> str:
    suffix = feature_key.split(".", 1)[1]
    if suffix == "lak_x":
        return "mgga_exc/mgga_x_lak.mpl"
    if suffix == "lak_c":
        return "mgga_exc/mgga_c_lak.mpl"
    return f"mgga_exc/mgga_xc_ml25_{suffix}.mpl"


def maple_binding_for_key(feature_key: str) -> str:
    suffix = feature_key.split(".", 1)[1]
    return f"{suffix}_h"


def parse_maple_formula_entries(path: Path) -> list[tuple[str, str]]:
    entries: list[tuple[str, str]] = []
    current_include: str | None = None
    for raw_line in read_text(path).splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("ifxc_ml25_features := ["):
            break
        include_match = INCLUDE_RE.match(line)
        if include_match:
            current_include = include_match.group(1)
            continue
        assign_match = ASSIGN_RE.match(line)
        if assign_match and current_include is not None:
            entries.append((current_include, assign_match.group(1)))
            current_include = None
    return entries


def parse_maple_keys(path: Path) -> list[tuple[int, str]]:
    keys: list[tuple[int, str]] = []
    in_vector = False
    for raw_line in read_text(path).splitlines():
        line = raw_line.strip()
        if line.startswith("ifxc_ml25_features := ["):
            in_vector = True
            continue
        if not in_vector:
            continue
        if line.startswith("]:"):
            break
        if not line or line.startswith("#"):
            continue
        match = MAPLE_KEY_RE.search(raw_line)
        if match:
            keys.append((int(match.group(1)), match.group(2)))
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
    repo_root = source_path.parent.parent

    source_text = read_text(source_path)
    if "ifxc_feature_set_key := \"ml25\":" not in source_text:
        return fail("Maple source does not advertise feature set ml25")
    if "ifxc_feature_count := 66:" not in source_text:
        return fail("Maple source does not advertise a 66-feature vector")

    feature_defs = parse_feature_def(feature_def_path)
    if len(feature_defs) != 66:
        return fail(f"expected 66 feature definitions, found {len(feature_defs)}")
    if [item[0] for item in feature_defs] != list(range(66)):
        return fail("feature definition indexes are not contiguous")

    maple_keys = parse_maple_keys(source_path)
    if len(maple_keys) != 66:
        return fail(f"expected 66 Maple feature entries, found {len(maple_keys)}")
    if [item[0] for item in maple_keys] != list(range(66)):
        return fail("Maple feature indexes are not contiguous")

    maple_entries = parse_maple_formula_entries(source_path)
    if len(maple_entries) != 66:
        return fail(f"expected 66 Maple feature includes, found {len(maple_entries)}")

    feature_keys = [item[2] for item in feature_defs]
    maple_feature_keys = [item[1] for item in maple_keys]
    if feature_keys != maple_feature_keys:
        return fail("feature def keys and Maple feature keys differ")

    expected_includes = [maple_include_for_key(key) for key in feature_keys]
    expected_bindings = [maple_binding_for_key(key) for key in feature_keys]
    actual_includes = [item[0] for item in maple_entries]
    actual_bindings = [item[1] for item in maple_entries]
    if actual_includes != expected_includes:
        return fail("Maple include list differs from feature definitions")
    if actual_bindings != expected_bindings:
        return fail("Maple binding list differs from feature definitions")

    for include_path in expected_includes:
        if not (source_path.parent / include_path).exists():
            return fail(f"missing Maple include: {include_path}")
    for support_path in REQUIRED_SUPPORT_FILES:
        if not (repo_root / support_path).exists():
            return fail(f"missing Maple support file: {support_path}")
    for generated_path in FORBIDDEN_GENERATED_METADATA:
        if (repo_root / generated_path).exists():
            return fail(f"redundant generated metadata file remains: {generated_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
