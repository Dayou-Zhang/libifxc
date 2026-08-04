#!/usr/bin/env python3
# SPDX-License-Identifier: MPL-2.0
"""IFXC vector-valued MGGA Maple2C generator scaffold."""

from __future__ import annotations

from maple2c_lib.mgga import work_mgga_exc


def work_if_mgga(params):
  """Generate an IFXC vector MGGA feature kernel.

  This is the Milestone-3 scaffold: it recognizes the IFXC Maple source type,
  validates IFXC-specific source headers, and emits into src/maple2c/if_mgga.
  The MGGA Maple2C backend already handles vector-valued feature derivatives.
  This wrapper validates IFXC-specific source headers and emits into
  src/maple2c/if_mgga.
  """

  if params["feature_set"] is None:
    raise SystemExit("if_mgga Maple sources require '(* feature_set: ... *)'")

  if params["n_features"] <= 0:
    raise SystemExit("if_mgga Maple sources require '(* nfeatures: N *)'")

  if params["source_max_order"] is not None and params["maxorder"] > params["source_max_order"]:
    raise SystemExit(
      "requested maxorder {} exceeds source max_order {}".format(
        params["maxorder"], params["source_max_order"]))

  if params["source_variables"] and params["source_variables"] != ["rho", "sigma", "tau"]:
    raise SystemExit("first-release if_mgga sources must declare variables: rho sigma tau")

  routed = params.copy()
  routed["functype"] = "if_mgga"
  routed["feature_batches"] = params["maxorder"] > 0
  routed["feature_batch_size"] = 4
  routed["max_workers"] = 1
  work_mgga_exc(routed)
