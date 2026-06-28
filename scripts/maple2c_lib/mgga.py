#!/usr/bin/env python3

# Copyright (C) 2021 M.A.L. Marques
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from maple2c_lib.utils import *

# these are the variables that the functional depends on
variables = ["rho_0_", "rho_1_", "sigma_0_", "sigma_1_", "sigma_2_", "lapl_0_", "lapl_1_", "tau_0_", "tau_1_"]
  
# get arguments of the functions
input_args  = "const double *rho, const double *sigma, const double *lapl, const double *tau"
output_args = "xc_mgga_out_params *out"

# the definition of the derivatives that libxc transmits to the calling program
partials = [
  ["zk"],
  ["vrho", "vsigma", "vlapl", "vtau"],
  ["v2rho2", "v2rhosigma", "v2rholapl", "v2rhotau", "v2sigma2",
   "v2sigmalapl", "v2sigmatau", "v2lapl2", "v2lapltau", "v2tau2"],
  ["v3rho3", "v3rho2sigma", "v3rho2lapl", "v3rho2tau", "v3rhosigma2",
   "v3rhosigmalapl", "v3rhosigmatau", "v3rholapl2", "v3rholapltau",
   "v3rhotau2", "v3sigma3", "v3sigma2lapl", "v3sigma2tau", "v3sigmalapl2",
   "v3sigmalapltau", "v3sigmatau2", "v3lapl3", "v3lapl2tau", "v3lapltau2",
   "v3tau3"],
  ["v4rho4", "v4rho3sigma", "v4rho3lapl", "v4rho3tau", "v4rho2sigma2",
   "v4rho2sigmalapl", "v4rho2sigmatau", "v4rho2lapl2", "v4rho2lapltau",
   "v4rho2tau2", "v4rhosigma3", "v4rhosigma2lapl", "v4rhosigma2tau",
   "v4rhosigmalapl2", "v4rhosigmalapltau", "v4rhosigmatau2",
   "v4rholapl3", "v4rholapl2tau", "v4rholapltau2", "v4rhotau3",
   "v4sigma4", "v4sigma3lapl", "v4sigma3tau", "v4sigma2lapl2",
   "v4sigma2lapltau", "v4sigma2tau2", "v4sigmalapl3", "v4sigmalapl2tau",
   "v4sigmalapltau2", "v4sigmatau3", "v4lapl4", "v4lapl3tau",
   "v4lapl2tau2", "v4lapltau3", "v4tau4",
  ]
]

#####################################################################
def expand_feature_derivatives(derivatives, n_features):
  expanded = []

  for der_order in derivatives:
    expanded.append([])
    for ifeat in range(n_features):
      for order, name in der_order:
        expanded[-1].append([order.copy(), flatten_feature_output_name(name, n_features, ifeat)])

  return expanded


def maple_define_feature_derivatives(variables, derivatives, func, n_features):
  out_derivatives = [[]]
  out_cgeneration = [[] for _ in derivatives]
  realvars = ", ".join(variables)

  for der_order in derivatives[1:]:
    out_derivatives.append([])
    for order, name in der_order:
      to_derive = order.copy()

      if all(v == 0 for v in to_derive):
        break

      for i_to_derive in range(len(to_derive)):
        if to_derive[i_to_derive] != 0:
          to_derive[i_to_derive] -= 1
          break

      vars = "v" + ", v".join(str(i) for i in range(len(order)))
      if all(v == 0 for v in to_derive):
        f_to_derive = "{}({})".format(func, vars)
      else:
        f_to_derive = "d{}d{}({})".format(func, "".join(str(i) for i in to_derive), vars)

      varname = "d{}d{}".format(func, "".join(str(i) for i in order))
      out_derivatives[-1].append("{} := ({}) -> [seq(eval(diff(op(i, {}), v{})), i = 1..{})]:\n\n".format(
        varname, vars, f_to_derive, i_to_derive, n_features)
      )

      for ifeat in range(n_features):
        out_cgeneration[len(out_derivatives) - 1].append(
          "{} = op({}, {}({}))".format(
            flatten_feature_output_name(name, n_features, ifeat), ifeat + 1, varname, realvars))

  return ["".join(der_order) for der_order in out_derivatives], out_cgeneration


def maple_define_feature_batch_derivatives(variables, batch_derivatives, func, n_features, maxorder):
  vars = "v" + ", v".join(str(i) for i in range(len(variables)))
  needed_orders = [[] for _ in range(maxorder + 1)]
  seen = set()
  pending = []

  for order, name in batch_derivatives:
    _, varorder = parse_output_name(name)
    pending.append((order.copy(), varorder % n_features))

  while len(pending) > 0:
    order, feature_index = pending.pop()
    order_key = (tuple(order), feature_index)
    if order_key in seen or all(v == 0 for v in order):
      continue

    seen.add(order_key)
    needed_orders[sum(order)].append((order, feature_index))

    previous = order.copy()
    for i_to_derive in range(len(previous)):
      if previous[i_to_derive] != 0:
        previous[i_to_derive] -= 1
        break

    if not all(v == 0 for v in previous):
      pending.append((previous, feature_index))

  out_derivatives = [""] * (maxorder + 1)
  for order_group in needed_orders:
    order_group.sort()

  for order_index in range(1, maxorder + 1):
    for order, feature_index in needed_orders[order_index]:
      previous = order.copy()
      for i_to_derive in range(len(previous)):
        if previous[i_to_derive] != 0:
          previous[i_to_derive] -= 1
          break

      if all(v == 0 for v in previous):
        f_to_derive = "op({}, {}({}))".format(feature_index + 1, func, vars)
      else:
        f_to_derive = "d{}d{}_f{}({})".format(
          func, "".join(str(i) for i in previous), feature_index, vars)

      varname = "d{}d{}_f{}".format(func, "".join(str(i) for i in order), feature_index)
      out_derivatives[order_index] += "{} := ({}) -> eval(diff({}, v{})):\n\n".format(
        varname, vars, f_to_derive, next(i for i, value in enumerate(order) if value != 0))

  return out_derivatives


def maple_define_feature_batches(variant_code, maple_code, outputs, derivatives, variables, func, n_features, feature_batch_size=0):
  batches = []

  def group_by_key(items, key_getter):
    groups = []
    current_group = []
    current_key = None

    for item in items:
      key = key_getter(item)
      if current_key is None or key == current_key:
        current_group.append(item)
      else:
        groups.append(current_group)
        current_group = [item]
      current_key = key

    if len(current_group) > 0:
      groups.append(current_group)

    return groups

  for order in range(len(derivatives)):
    if len(outputs[order]) == 0:
      continue

    derivative_lookup = {
      name: der_order.copy() for der_order, name in derivatives[order]
    }
    output_groups = group_by_key(
      outputs[order],
      lambda item: tuple(derivative_lookup[item.split(" =", 1)[0].strip()]),
    )

    for output_group in output_groups:
      if feature_batch_size > 0 and order > 0:
        output_batches = [
          output_group[i:i + feature_batch_size]
          for i in range(0, len(output_group), feature_batch_size)
        ]
      else:
        output_batches = [output_group]

      for batch_outputs in output_batches:
        batch_derivatives = []
        for output in batch_outputs:
          output_name = output.split(" =", 1)[0].strip()
          batch_derivatives.append([derivative_lookup[output_name].copy(), output_name])

        batch_defs = maple_define_feature_batch_derivatives(
          variables, batch_derivatives, func, n_features, len(derivatives) - 1)
        batch_c_outputs = batch_outputs
        if order > 0:
          realvars = ", ".join(variables)
          batch_c_outputs = []
          for der_order, output_name in batch_derivatives:
            _, varorder = parse_output_name(output_name)
            feature_index = varorder % n_features
            batch_c_outputs.append("{} = d{}d{}_f{}({})".format(
              output_name, func, "".join(str(i) for i in der_order), feature_index, realvars))

        batches.append({
          "start_order": order,
          "label": "order {} {}".format(
            order, batch_outputs[0].split(" =", 1)[0].strip()),
          "derivatives": [batch_derivatives],
          "code": '''
{}

{}
{}
C([{}], optimize, deducetypes=false):

'''.format(
          variant_code,
          "".join(batch_defs[1:]),
          maple_code,
          ", ".join(batch_c_outputs),
        )
        })

  return {"batches": batches}


def maple_define_feature_roots(params, variables):
  maple_zk = []
  realvars = ", ".join(variables)
  simplify_line = ""

  if params["simplify_begin"] != "":
    simplify_line = "  feat := [seq(simplify(op(i, feat), symbolic), i = 1..{})]:\n".format(params["n_features"])

  maple_roots = '''
mzkf := proc(r0, r1, s0, s1, s2, l0, l1, tau0, tau1)
  local feat, i:

  feat := f(r_ws(dens(r0, r1)), zeta(r0, r1), xt(r0, r1, s0, s1, s2), xs0(r0, r1, s0, s2), xs1(r0, r1, s0, s2), u0(r0, r1, l0, l1), u1(r0, r1, l0, l1), t0(r0, r1, tau0, tau1), t1(r0, r1, tau0, tau1)):
{}  return feat:
end proc:

mff := proc(r0, r1, s0, s1, s2, l0, l1, tau0, tau1)
  local feat, i:

  feat := mzkf(r0, r1, s0, s1, s2, l0, l1, tau0, tau1):
  return [seq(dens(r0, r1)*op(i, feat), i = 1..{})]:
end proc:

'''.format(simplify_line, params["n_features"])

  for ifeat in range(params["n_features"]):
    maple_zk.append(" zk_{}_ = op({}, mzkf({}))".format(ifeat, ifeat + 1, realvars))

  return maple_roots, maple_zk


#####################################################################
def work_mgga_exc(params):
  '''Process a MGGA functional for the energy'''

  derivatives = partials_to_derivatives(params, "mgga", partials)
  if params.get("source_variables") == ["rho", "sigma", "tau"]:
    derivatives = [
      [der for der in der_order if "lapl" not in der[1]]
      for der_order in derivatives
    ]

  if params["n_features"] > 0:
    base_derivatives_unpol = filter_unpolarized_derivatives(derivatives)
    derivatives_unpol = expand_feature_derivatives(base_derivatives_unpol, params["n_features"])
    derivatives_pol = expand_feature_derivatives(derivatives, params["n_features"])

    der_def_unpol, out_c_unpol = maple_define_feature_derivatives(
      variables, base_derivatives_unpol, "mff", params["n_features"])
    der_def_pol, out_c_pol = maple_define_feature_derivatives(
      variables, derivatives, "mff", params["n_features"])

    maple_roots, maple_zk = maple_define_feature_roots(params, variables)

    maple_code  = '''
# zk is energy per unit particle
{}
$include <util.mpl>
'''.format(maple_roots)

    if params["feature_batches"]:
      outputs_unpol = [maple_zk] + out_c_unpol[1:]
      outputs_pol = [maple_zk] + out_c_pol[1:]

      variants = {
        "unpol": maple_define_feature_batches(
          '''
dens := (r0, r1) -> r0:
zeta := (r0, r1) -> 0:
xs0  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0/4)/((r0/2)^(1 + 1/DIMENSIONS)):
xs1  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0/4)/((r0/2)^(1 + 1/DIMENSIONS)):
xt   := (r0, r1, sigma0, sigma1, sigma2) -> sqrt(sigma0)/r0^(1 + 1/DIMENSIONS):
u0   := (r0, r1, l0, l1) -> (l0/2)/((r0/2)^(1 + 2/DIMENSIONS)):
u1   := (r0, r1, l0, l1) -> (l0/2)/((r0/2)^(1 + 2/DIMENSIONS)):
t0   := (r0, r1, tau0, tau1) -> (tau0/2)/((r0/2)^(1 + 2/DIMENSIONS)):
t1   := (r0, r1, tau0, tau1) -> (tau0/2)/((r0/2)^(1 + 2/DIMENSIONS)):
''',
          maple_code,
          outputs_unpol,
          derivatives_unpol,
          variables,
          "mff",
          params["n_features"],
          params.get("feature_batch_size", 0)
        ),

        "pol": maple_define_feature_batches(
          '''
dens := (r0, r1) -> r0 + r1:
zeta := (r0, r1) -> (r0 - r1)/(r0 + r1):
xs0  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0)/r0^(1 + 1/DIMENSIONS):
xs1  := (r0, r1, sigma0, sigma2) -> sqrt(sigma2)/r1^(1 + 1/DIMENSIONS):
xt   := (r0, r1, sigma0, sigma1, sigma2) -> sqrt(sigma0 + 2*sigma1 + sigma2)/(r0 + r1)^(1 + 1/DIMENSIONS):
u0   := (r0, r1, l0, l1) -> l0/(r0^(1 + 2/DIMENSIONS)):
u1   := (r0, r1, l0, l1) -> l1/(r1^(1 + 2/DIMENSIONS)):
t0   := (r0, r1, tau0, tau1) -> tau0/(r0^(1 + 2/DIMENSIONS)):
t1   := (r0, r1, tau0, tau1) -> tau1/(r1^(1 + 2/DIMENSIONS)):
''',
          maple_code,
          outputs_pol,
          derivatives_pol,
          variables,
          "mff",
          params["n_features"],
          params.get("feature_batch_size", 0)
        )
      }
    else:
      der_def_unpol = "".join(der_def_unpol)
      der_def_pol = "".join(der_def_pol)

      maple_zk = ", ".join(maple_zk)
      out_c_unpol = ", ".join([line for lines in out_c_unpol for line in lines])
      out_c_pol = ", ".join([line for lines in out_c_pol for line in lines])

      if out_c_unpol != "":
        out_c_unpol = ", " + out_c_unpol
      if out_c_pol != "":
        out_c_pol = ", " + out_c_pol

      variants = {
        "unpol": {
          "code": '''
dens := (r0, r1) -> r0:
zeta := (r0, r1) -> 0:
xs0  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0/4)/((r0/2)^(1 + 1/DIMENSIONS)):
xs1  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0/4)/((r0/2)^(1 + 1/DIMENSIONS)):
xt   := (r0, r1, sigma0, sigma1, sigma2) -> sqrt(sigma0)/r0^(1 + 1/DIMENSIONS):
u0   := (r0, r1, l0, l1) -> (l0/2)/((r0/2)^(1 + 2/DIMENSIONS)):
u1   := (r0, r1, l0, l1) -> (l0/2)/((r0/2)^(1 + 2/DIMENSIONS)):
t0   := (r0, r1, tau0, tau1) -> (tau0/2)/((r0/2)^(1 + 2/DIMENSIONS)):
t1   := (r0, r1, tau0, tau1) -> (tau0/2)/((r0/2)^(1 + 2/DIMENSIONS)):

{}
{}
C([{}{}], optimize, deducetypes=false):

'''.format(der_def_unpol, maple_code, maple_zk, out_c_unpol),
          "derivatives": derivatives_unpol,
        },

        "pol": {
          "code": '''
dens := (r0, r1) -> r0 + r1:
zeta := (r0, r1) -> (r0 - r1)/(r0 + r1):
xs0  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0)/r0^(1 + 1/DIMENSIONS):
xs1  := (r0, r1, sigma0, sigma2) -> sqrt(sigma2)/r1^(1 + 1/DIMENSIONS):
xt   := (r0, r1, sigma0, sigma1, sigma2) -> sqrt(sigma0 + 2*sigma1 + sigma2)/(r0 + r1)^(1 + 1/DIMENSIONS):
u0   := (r0, r1, l0, l1) -> l0/(r0^(1 + 2/DIMENSIONS)):
u1   := (r0, r1, l0, l1) -> l1/(r1^(1 + 2/DIMENSIONS)):
t0   := (r0, r1, tau0, tau1) -> tau0/(r0^(1 + 2/DIMENSIONS)):
t1   := (r0, r1, tau0, tau1) -> tau1/(r1^(1 + 2/DIMENSIONS)):

{}
{}
C([{}{}], optimize, deducetypes=false):

'''.format(der_def_pol, maple_code, maple_zk, out_c_pol),
          "derivatives": derivatives_pol,
        }
      }
  else:
    der_def, out_c = maple_define_derivatives(variables, derivatives, "mf")

    out_c = ", ".join(out_c)
    if out_c != "": out_c = ", " + out_c

    # we join all the pieces
    maple_code  = '''
# zk is energy per unit particle
mzk  := (r0, r1, s0, s1, s2, l0, l1, tau0, tau1) -> \\
  {} + \\
    f(r_ws(dens(r0, r1)), zeta(r0, r1), xt(r0, r1, s0, s1, s2), xs0(r0, r1, s0, s2), xs1(r0, r1, s0, s2), u0(r0, r1, l0, l1), u1(r0, r1, l0, l1), t0(r0, r1, tau0, tau1), t1(r0, r1, tau0, tau1)) \\
  {} :

  (* mf is energy per unit volume *)
  mf   := (r0, r1, s0, s1, s2, l0, l1, tau0, tau1) -> eval(dens(r0, r1)*mzk(r0, r1, s0, s1, s2, l0, l1, tau0, tau1)):

$include <util.mpl>
'''.format(params["simplify_begin"], params["simplify_end"])

    maple_zk = " zk_0_ = mzk(" + ", ".join(variables) + ")"

    # we build 2 variants of the functional, for unpolarized, and polarized densities
    variants = {
      "unpol": '''
dens := (r0, r1) -> r0:
zeta := (r0, r1) -> 0:
xs0  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0/4)/((r0/2)^(1 + 1/DIMENSIONS)):
xs1  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0/4)/((r0/2)^(1 + 1/DIMENSIONS)):
xt   := (r0, r1, sigma0, sigma1, sigma2) -> sqrt(sigma0)/r0^(1 + 1/DIMENSIONS):
u0   := (r0, r1, l0, l1) -> (l0/2)/((r0/2)^(1 + 2/DIMENSIONS)):
u1   := (r0, r1, l0, l1) -> (l0/2)/((r0/2)^(1 + 2/DIMENSIONS)):
t0   := (r0, r1, tau0, tau1) -> (tau0/2)/((r0/2)^(1 + 2/DIMENSIONS)):
t1   := (r0, r1, tau0, tau1) -> (tau0/2)/((r0/2)^(1 + 2/DIMENSIONS)):

{}

{}
C([{}{}], optimize, deducetypes=false):

'''.format(der_def, maple_code, maple_zk, out_c),

      "pol": '''
dens := (r0, r1) -> r0 + r1:
zeta := (r0, r1) -> (r0 - r1)/(r0 + r1):
xs0  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0)/r0^(1 + 1/DIMENSIONS):
xs1  := (r0, r1, sigma0, sigma2) -> sqrt(sigma2)/r1^(1 + 1/DIMENSIONS):
xt   := (r0, r1, sigma0, sigma1, sigma2) -> sqrt(sigma0 + 2*sigma1 + sigma2)/(r0 + r1)^(1 + 1/DIMENSIONS):
u0   := (r0, r1, l0, l1) -> l0/(r0^(1 + 2/DIMENSIONS)):
u1   := (r0, r1, l0, l1) -> l1/(r1^(1 + 2/DIMENSIONS)):
t0   := (r0, r1, tau0, tau1) -> tau0/(r0^(1 + 2/DIMENSIONS)):
t1   := (r0, r1, tau0, tau1) -> tau1/(r1^(1 + 2/DIMENSIONS)):

{}

{}
C([{}{}], optimize, deducetypes=false):

'''.format(der_def, maple_code, maple_zk, out_c)
    }

  maple2c_run(params, variables, derivatives, variants, 0, input_args, output_args)


#####################################################################
def work_mgga_vxc(params):
  '''Process a MGGA functional for the potential'''

  all_derivatives = partials_to_derivatives(params, "mgga", partials)

  derivatives, derivatives1, derivatives2 = filter_vxc_derivatives(all_derivatives)
  
  # we obtain the missing pieces for maple
  # unpolarized calculation
  der_def_unpol, out_c_unpol = maple_define_derivatives(variables, derivatives1, "mf0")
  out_c_unpol = ", ".join(out_c_unpol)
  if out_c_unpol != "": out_c_unpol = ", " + out_c_unpol

  # polarized calculation
  der_def_pol1, out_c_pol1 = maple_define_derivatives(variables, derivatives1, "mf0")
  der_def_pol2, out_c_pol2 = maple_define_derivatives(variables, derivatives2, "mf1")
  
  der_def_pol = der_def_pol1 + der_def_pol2
  out_c_pol   = ", ".join(sorted(out_c_pol1 + out_c_pol2, key=sort_alphanumerically))
  if out_c_pol != "": out_c_pol = ", " + out_c_pol
  
  # we join all the pieces
  maple_code  = '''
mzk  := (r0, r1, s0, s1, s2, l0, l1, tau0, tau1) -> \\
  {} + \\
    f(r_ws(dens(r0, r1)), zeta(r0, r1), xt(r0, r1, s0, s1, s2), xs0(r0, r1, s0, s2), xs1(r0, r1, s0, s2), u0(r0, r1, l0, l1), u1(r0, r1, l0, l1), t0(r0, r1, tau0, tau1), t1(r0, r1, tau0, tau1)) \\
  {} :

(* mf is the up potential *)
mf0   := (r0, r1, s0, s1, s2, l0, l1, tau0, tau1) -> eval(mzk(r0, r1, s0, s1, s2, l0, l1, tau0, tau1)):
mf1   := (r0, r1, s0, s1, s2, l0, l1, tau0, tau1) -> eval(mzk(r1, r0, s2, s1, s0, l1, l0, tau1, tau0)):

$include <util.mpl>
'''.format(params["simplify_begin"], params["simplify_end"])
  
  maple_vrho0 = " vrho_0_ = mf0(" + ", ".join(variables) + ")"
  maple_vrho1 = " vrho_1_ = mf1(" + ", ".join(variables) + ")"

  # we build 2 variants of the functional, for unpolarized, and polarized densities
  variants = {
    "unpol": '''
dens := (r0, r1) -> r0:
zeta := (r0, r1) -> 0:
xs0  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0/4)/((r0/2)^(1 + 1/DIMENSIONS)):
xs1  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0/4)/((r0/2)^(1 + 1/DIMENSIONS)):
xt   := (r0, r1, sigma0, sigma1, sigma2) -> sqrt(sigma0)/r0^(1 + 1/DIMENSIONS):
u0   := (r0, r1, l0, l1) -> (l0/2)/((r0/2)^(1 + 2/DIMENSIONS)):
u1   := (r0, r1, l0, l1) -> (l0/2)/((r0/2)^(1 + 2/DIMENSIONS)):
t0   := (r0, r1, tau0, tau1) -> (tau0/2)/((r0/2)^(1 + 2/DIMENSIONS)):
t1   := (r0, r1, tau0, tau1) -> (tau0/2)/((r0/2)^(1 + 2/DIMENSIONS)):

{}

{}
C([{}{}], optimize, deducetypes=false):

'''.format(der_def_unpol, maple_code, maple_vrho0, out_c_unpol),

    "pol": '''
dens := (r0, r1) -> r0 + r1:
zeta := (r0, r1) -> (r0 - r1)/(r0 + r1):
xs0  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0)/r0^(1 + 1/DIMENSIONS):
xs1  := (r0, r1, sigma0, sigma2) -> sqrt(sigma2)/r1^(1 + 1/DIMENSIONS):
xt   := (r0, r1, sigma0, sigma1, sigma2) -> sqrt(sigma0 + 2*sigma1 + sigma2)/(r0 + r1)^(1 + 1/DIMENSIONS):
u0   := (r0, r1, l0, l1) -> l0/(r0^(1 + 2/DIMENSIONS)):
u1   := (r0, r1, l0, l1) -> l1/(r1^(1 + 2/DIMENSIONS)):
t0   := (r0, r1, tau0, tau1) -> tau0/(r0^(1 + 2/DIMENSIONS)):
t1   := (r0, r1, tau0, tau1) -> tau1/(r1^(1 + 2/DIMENSIONS)):

{}

{}
C([{}, {}{}], optimize, deducetypes=false):

'''.format(der_def_pol, maple_code, maple_vrho0, maple_vrho1, out_c_pol)
  }

  maple2c_run(params, variables, derivatives, variants, 1, input_args, output_args)
