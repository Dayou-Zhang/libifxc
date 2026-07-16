#!/usr/bin/env python3
"""Validate ML26 manifest metadata without regenerating formula code."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


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
    source = source_path.read_text(encoding="utf-8")
    feature_def = feature_def_path.read_text(encoding="utf-8")

    for directive in (
        "(* type: if_mgga *)",
        "(* feature_set: ml26 *)",
        "(* nfeatures: 69 *)",
        "(* max_order: 0 *)",
        "(* variables: rho sigma tau *)",
        '$include "mgga_xc_ml25.mpl"',
        'ifxc_feature_set_key := "ml26":',
        "ifxc_feature_count := 69:",
    ):
        if directive not in source:
            return fail(f"missing ML26 Maple directive: {directive}")

    labels = (
        "ml26.cs1_same_spin_constant",
        "ml26.cs1_same_spin_gradient",
        "ml26.cs1_opposite_spin",
    )
    for index, label in enumerate(labels, start=66):
        if f'"{label}"' not in source:
            return fail(f"missing ML26 Maple feature label: {label}")
        if f"IFXC_ML26_FEATURE({index}," not in feature_def:
            return fail(f"missing ML26 feature definition index {index}")
        if f'"{label}"' not in feature_def:
            return fail(f"ML26 feature definition label differs: {label}")

    if '#include "ifxc_ml25_features.def"' not in feature_def:
        return fail("ML26 feature metadata does not reuse the ML25 prefix")

    generated_path = repo_root / "src/maple2c/if_mgga/mgga_xc_ml26.c"
    if not generated_path.exists():
        return fail("missing generated ML26 Maple2C output")
    generated = generated_path.read_text(encoding="utf-8")
    for marker in (
        "Type of functional: if_mgga",
        "#define ifxc_maple2c_order 0",
        "out->zk[ip*p->dim.zk + 68]",
    ):
        if marker not in generated:
            return fail(f"generated ML26 output is missing marker: {marker}")
    if "Error," in generated:
        return fail("Maple error text remains in generated ML26 output")
    if "out->vrho" in generated or "out->vsigma" in generated or "out->vtau" in generated:
        return fail("order-0 ML26 output unexpectedly contains derivatives")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
