#!/usr/bin/env python3

# Copyright (C) 2021 M.A.L. Marques
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import sys, os, re, subprocess
from concurrent.futures import ThreadPoolExecutor

# we need this in a couple of places
der_name = ("EXC", "VXC", "FXC", "KXC", "LXC", "MXC")

#####################################################################
# sort by character and then by number
def sort_alphanumerically(x):
  res = re.match(r"([^_]+)_([0-9]+)_", x)
  return (res.group(1), int(res.group(2)))


def parse_output_name(name):
  res = re.match(r"([^_]+)_([0-9]+)_", name)
  if res is None:
    raise ValueError("Could not parse output name '{}'".format(name))

  return res.group(1), int(res.group(2))


def flatten_feature_output_name(name, n_features, feature_number):
  varname, varorder = parse_output_name(name)
  return "{}_{}_".format(varname, varorder*n_features + feature_number)


#####################################################################
def partials_to_derivatives(params, func_type, partials):
  derivatives = []
  for order in range(params["maxorder"] + 1):
    derivatives.append([])
    for der in partials[order]:
      derivatives[order].extend(enumerate_spin_partials(der, func_type))

  return derivatives


#####################################################################
def filter_vxc_derivatives(all_derivatives):
  '''This separates the derivatives of vxc into derivatives of vxc_0
  and vxc_1. All other derivatives (e.g. vsigma) are ignored.'''

  derivatives  = []
  derivatives1 = []
  derivatives2 = []

  for order in range(len(all_derivatives) - 1):
    derivatives.append([])
    derivatives1.append([])
    derivatives2.append([])
    
    for der in all_derivatives[order + 1]:
      if der[0][0] > 0:
        der[0][0] -= 1
        derivatives1[order].append(der)
        derivatives [order].append(der)
      elif der[0][1] > 0:
        der[0][1] -= 1
        derivatives2[order].append(der)
        derivatives [order].append(der)

  return derivatives, derivatives1, derivatives2


def filter_unpolarized_derivatives(all_derivatives):
  derivatives = []

  for der_order in all_derivatives:
    derivatives.append([])
    for der in der_order:
      if parse_output_name(der[1])[1] == 0:
        derivatives[-1].append([der[0].copy(), der[1]])

  return derivatives


def enumerate_spin_partials(derivative, func_type):
  '''Given the name of a derivative (such as 'v2rho2')
  and a functional type ("lda", "gga", "mgga"), return all spin 
  variants such as

  [
    [[2, 0], 'v2rho2_0_'], 
    [[1, 1], 'v2rho2_1_'], 
    [[0, 2], 'v2rho2_2_']
  ]
  '''

  words = ("rho", "sigma", "lapl", "tau")
  partials = {
    "rho"   :  [[[0, 0, 0, 0, 0, 0, 0, 0, 0]],    # 0th-order
                [[1, 0, 0, 0, 0, 0, 0, 0, 0],     # 1st-order
                 [0, 1, 0, 0, 0, 0, 0, 0, 0]],
                [[2, 0, 0, 0, 0, 0, 0, 0, 0],     # 2nd-order
                 [1, 1, 0, 0, 0, 0, 0, 0, 0],
                 [0, 2, 0, 0, 0, 0, 0, 0, 0]],
                [[3, 0, 0, 0, 0, 0, 0, 0, 0],     # 3rd-order
                 [2, 1, 0, 0, 0, 0, 0, 0, 0],
                 [1, 2, 0, 0, 0, 0, 0, 0, 0],
                 [0, 3, 0, 0, 0, 0, 0, 0, 0]],
                [[4, 0, 0, 0, 0, 0, 0, 0, 0],     # 4th-order
                 [3, 1, 0, 0, 0, 0, 0, 0, 0],
                 [2, 2, 0, 0, 0, 0, 0, 0, 0],
                 [1, 3, 0, 0, 0, 0, 0, 0, 0],
                 [0, 4, 0, 0, 0, 0, 0, 0, 0]],
    ],
    "sigma" :  [[[0, 0, 0, 0, 0, 0, 0, 0, 0]],    # 0th-order
                [[0, 0, 1, 0, 0, 0, 0, 0, 0],     # 1st-order
                 [0, 0, 0, 1, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 1, 0, 0, 0, 0]],
                [[0, 0, 2, 0, 0, 0, 0, 0, 0],     # 2nd-order
                 [0, 0, 1, 1, 0, 0, 0, 0, 0],
                 [0, 0, 1, 0, 1, 0, 0, 0, 0],
                 [0, 0, 0, 2, 0, 0, 0, 0, 0],
                 [0, 0, 0, 1, 1, 0, 0, 0, 0],
                 [0, 0, 0, 0, 2, 0, 0, 0, 0]],
                [[0, 0, 3, 0, 0, 0, 0, 0, 0],     # 3rd-order
                 [0, 0, 2, 1, 0, 0, 0, 0, 0],
                 [0, 0, 2, 0, 1, 0, 0, 0, 0],
                 [0, 0, 1, 2, 0, 0, 0, 0, 0],
                 [0, 0, 1, 1, 1, 0, 0, 0, 0],
                 [0, 0, 1, 0, 2, 0, 0, 0, 0],
                 [0, 0, 0, 3, 0, 0, 0, 0, 0],
                 [0, 0, 0, 2, 1, 0, 0, 0, 0],
                 [0, 0, 0, 1, 2, 0, 0, 0, 0],
                 [0, 0, 0, 0, 3, 0, 0, 0, 0]],
                [[0, 0, 4, 0, 0, 0, 0, 0, 0],     # 4th-order
                 [0, 0, 3, 1, 0, 0, 0, 0, 0],
                 [0, 0, 3, 0, 1, 0, 0, 0, 0],
                 [0, 0, 2, 2, 0, 0, 0, 0, 0],
                 [0, 0, 2, 1, 1, 0, 0, 0, 0],
                 [0, 0, 2, 0, 2, 0, 0, 0, 0],
                 [0, 0, 1, 3, 0, 0, 0, 0, 0],
                 [0, 0, 1, 2, 1, 0, 0, 0, 0],
                 [0, 0, 1, 1, 2, 0, 0, 0, 0],
                 [0, 0, 1, 0, 3, 0, 0, 0, 0],
                 [0, 0, 0, 4, 0, 0, 0, 0, 0],
                 [0, 0, 0, 3, 1, 0, 0, 0, 0],
                 [0, 0, 0, 2, 2, 0, 0, 0, 0],
                 [0, 0, 0, 1, 3, 0, 0, 0, 0],
                 [0, 0, 0, 0, 4, 0, 0, 0, 0]]
    ],
    "lapl" :   [[[0, 0, 0, 0, 0, 0, 0, 0, 0]],    # 0th-order
                [[0, 0, 0, 0, 0, 1, 0, 0, 0],     # 1st-order
                 [0, 0, 0, 0, 0, 0, 1, 0, 0]],
                [[0, 0, 0, 0, 0, 2, 0, 0, 0],     # 2nd-order
                 [0, 0, 0, 0, 0, 1, 1, 0, 0],
                 [0, 0, 0, 0, 0, 0, 2, 0, 0]],
                [[0, 0, 0, 0, 0, 3, 0, 0, 0],     # 3rd-order
                 [0, 0, 0, 0, 0, 2, 1, 0, 0],
                 [0, 0, 0, 0, 0, 1, 2, 0, 0],
                 [0, 0, 0, 0, 0, 0, 3, 0, 0]],
                [[0, 0, 0, 0, 0, 4, 0, 0, 0],     # 4th-order
                 [0, 0, 0, 0, 0, 3, 1, 0, 0],
                 [0, 0, 0, 0, 0, 2, 2, 0, 0],
                 [0, 0, 0, 0, 0, 1, 3, 0, 0],
                 [0, 0, 0, 0, 0, 0, 4, 0, 0]]
    ],
    "tau"  :   [[[0, 0, 0, 0, 0, 0, 0, 0, 0]],    # 0th-order
                [[0, 0, 0, 0, 0, 0, 0, 1, 0],     # 1st-order
                 [0, 0, 0, 0, 0, 0, 0, 0, 1]],
                [[0, 0, 0, 0, 0, 0, 0, 2, 0],     # 2nd-order
                 [0, 0, 0, 0, 0, 0, 0, 1, 1],
                 [0, 0, 0, 0, 0, 0, 0, 0, 2]],
                [[0, 0, 0, 0, 0, 0, 0, 3, 0],     # 3rd-order
                 [0, 0, 0, 0, 0, 0, 0, 2, 1],
                 [0, 0, 0, 0, 0, 0, 0, 1, 2],
                 [0, 0, 0, 0, 0, 0, 0, 0, 3]],
                [[0, 0, 0, 0, 0, 0, 0, 4, 0],     # 4th-order
                 [0, 0, 0, 0, 0, 0, 0, 3, 1],
                 [0, 0, 0, 0, 0, 0, 0, 2, 2],
                 [0, 0, 0, 0, 0, 0, 0, 1, 3],
                 [0, 0, 0, 0, 0, 0, 0, 0, 4]]
                ]
  }

  # finds out the order of each partial
  order = {}
  for word in words:
    order[word] = 0
  
    m = re.match(r'.*' + word + r'([0-9]*)', derivative)
    if m is not None:
      order[word] = 1 if m.group(1) == "" else int(m.group(1))

  # and this is the order of the derivative
  total_order = sum(order.values())

  max_n = {"lda": 2, "gga": 5, "mgga": 9}[func_type]

  all_derivatives = []
  der_n = 0
  for n_rho, p_rho in enumerate(partials["rho"][order["rho"]]):
    for n_sigma, p_sigma in enumerate(partials["sigma"][order["sigma"]]):
      for n_lapl, p_lapl in enumerate(partials["lapl"][order["lapl"]]):
        for n_tau, p_tau in enumerate(partials["tau"][order["tau"]]):
          # sum orders in all variables
            
          final_der = [0] * max_n
          for i in range(max_n):
            final_der[i] += \
                partials["rho"][order["rho"]][n_rho][i] + \
                partials["sigma"][order["sigma"]][n_sigma][i] + \
                partials["lapl"][order["lapl"]][n_lapl][i] + \
                partials["tau"][order["tau"]][n_tau][i]

          all_derivatives.append([final_der, derivative + "_" + str(der_n) + "_"])
          der_n += 1

  return all_derivatives


def maple_define_derivatives(variables, derivatives, func, name_transform=None, base_expr_builder=None):
  '''Generates maple code to define derivatives. Output is

  out_derivatives = "
dmfd10 := (v0, v1) ->  eval(diff(mf(v0, v1), v0)):

dmfd01 := (v0, v1) ->  eval(diff(mf(v0, v1), v1)):
...
"
  out_cgeneration = "[
    'vrho_0_ = dmfd10(rho_0_, rho_1_)', 
    'vrho_1_ = dmfd01(rho_0_, rho_1_)'
    ...
  ]
  
  '''
  
  out_derivatives = ""
  out_cgeneration = []

  realvars = ", ".join(variables)

  for der_order in derivatives:
    for der in der_order:
      order, name = der
      to_derive = order.copy()
      
      # is there something to do?
      if all([v == 0 for v in to_derive]):
        break

      # derivates are always defined as a first derivative of another derivative
      # so we have to find out what is this previous function to derive
      for i_to_derive in range(len(to_derive)):
        if to_derive[i_to_derive] != 0:
          to_derive[i_to_derive] -= 1
          break

      vars     = "v" + ", v".join(str(i) for i in range(len(order)))

      if all([v == 0 for v in to_derive]):
        if base_expr_builder is None:
          f_to_derive = func + "(" + vars + ")"
        else:
          f_to_derive = base_expr_builder(vars)
      else:
        f_to_derive = "d" + func + "d" + "".join(str(i) for i in to_derive) + "(" + vars + ")"

      # we build the expression of the derivative
      varname  = "d" + func + "d" + "".join(str(i) for i in order)
      derorder = ", ".join(str(i) for i in order)

      out_derivatives += varname + " := (" + vars + ") ->  eval(diff(" + \
        f_to_derive + ", v" + str(i_to_derive)  + ")):\n\n"

      out_name = name if name_transform is None else name_transform(name)
      out_cgeneration.append(out_name + " = " + varname + "(" + realvars + ")")

  return out_derivatives, out_cgeneration


def print_c_header(params, out):
  # Check for license expiry
  cmd = "echo -e 'quit;' | maple 2>&1 | grep License | head -n 1"
  license_expires = subprocess.check_output(cmd, shell=True).strip()
  if len(license_expires):
    print(f'Warning: {license_expires.decode("UTF-8")}')

  cmd = "echo -e 'quit;' | maple 2>&1 | grep Maple | head -n 1 | sed 's/^.*Maple/Maple/'"
  maple_version = subprocess.check_output(cmd, shell=True).strip()

  out.write('''/*
  This file was generated automatically with {}.
  Do not edit this file directly as it can be overwritten!!

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.

  Maple version     : {}
  Maple source      : {}
  Type of functional: {}
*/

#define ifxc_maple2c_order {}
'''.format(sys.argv[0], maple_version.decode(), params['maple_file'], params['functype'], params['maxorder']))


def strip_outer_parentheses(text):
  stripped = text.strip()

  while stripped.startswith("(") and stripped.endswith(")"):
    depth = 0
    balanced = True
    for i, char in enumerate(stripped):
      if char == "(":
        depth += 1
      elif char == ")":
        depth -= 1

      if depth == 0 and i != len(stripped) - 1:
        balanced = False
        break

    if not balanced:
      break

    stripped = stripped[1:-1].strip()

  return stripped


def split_function_args(text):
  depth = 0

  for i, char in enumerate(text):
    if char == "(":
      depth += 1
    elif char == ")":
      depth -= 1
    elif char == "," and depth == 0:
      return text[:i], text[i + 1:]

  return None


def replace_pow_calls(text):
  pow_replacements = {
    "0.1e1": lambda base: "({})".format(base),
    "-0.1e1": lambda base: "0.1e1/({})".format(base),
    "0.15e1": lambda base: "POW_3_2({})".format(base),
    "-0.15e1": lambda base: "0.1e1/POW_3_2({})".format(base),
    "0.2e1": lambda base: "POW_2({})".format(base),
    "-0.2e1": lambda base: "0.1e1/POW_2({})".format(base),
    "0.3e1": lambda base: "POW_3({})".format(base),
    "-0.3e1": lambda base: "0.1e1/POW_3({})".format(base),
    "0.5e0": lambda base: "sqrt({})".format(base),
    "-0.5e0": lambda base: "0.1e1/sqrt({})".format(base),
    "0.1e1/0.3e1": lambda base: "POW_1_3({})".format(base),
    "-0.1e1/0.3e1": lambda base: "0.1e1 / POW_1_3({})".format(base),
    "0.2e1/0.3e1": lambda base: "POW_2_3({})".format(base),
    "-0.2e1/0.3e1": lambda base: "0.1e1 / POW_2_3({})".format(base),
    "0.1e1/0.4e1": lambda base: "POW_1_4({})".format(base),
    "-0.1e1/0.4e1": lambda base: "0.1e1/POW_1_4({})".format(base),
    "0.333333333333333333.e0": lambda base: "POW_1_3({})".format(base),
    "-0.333333333333333333.e0": lambda base: "0.1e1 / POW_1_3({})".format(base),
    "0.666666666666666666.e0": lambda base: "POW_2_3({})".format(base),
    "-0.666666666666666666.e0": lambda base: "0.1e1 / POW_2_3({})".format(base),
    "0.1333333333333333333.e1": lambda base: "POW_4_3({})".format(base),
    "-0.1333333333333333333.e1": lambda base: "0.1e1 / POW_4_3({})".format(base),
    "0.4e1/0.3e1": lambda base: "POW_4_3({})".format(base),
    "-0.4e1/0.3e1": lambda base: "0.1e1 / POW_4_3({})".format(base),
    "0.1666666666666666666.e1": lambda base: "POW_5_3({})".format(base),
    "-0.1666666666666666666.e1": lambda base: "0.1e1 / POW_5_3({})".format(base),
    "0.5e1/0.3e1": lambda base: "POW_5_3({})".format(base),
    "-0.5e1/0.3e1": lambda base: "0.1e1 / POW_5_3({})".format(base),
    "0.2333333333333333333.e1": lambda base: "POW_7_3({})".format(base),
    "-0.2333333333333333333.e1": lambda base: "0.1e1 / POW_7_3({})".format(base),
    "0.7e1/0.3e1": lambda base: "POW_7_3({})".format(base),
    "-0.7e1/0.3e1": lambda base: "0.1e1 / POW_7_3({})".format(base),
  }

  def find_pow_start(source, start_index):
    while True:
      position = source.find("pow(", start_index)
      if position < 0:
        return -1

      if position == 0 or not (source[position - 1].isalnum() or source[position - 1] == "_"):
        return position

      start_index = position + 1

  def find_matching_paren(source, open_index):
    depth = 0
    for i in range(open_index, len(source)):
      if source[i] == "(":
        depth += 1
      elif source[i] == ")":
        depth -= 1
        if depth == 0:
          return i

    return -1

  result = []
  index = 0

  while index < len(text):
    start = find_pow_start(text, index)
    if start < 0:
      result.append(text[index:])
      break

    result.append(text[index:start])
    open_index = start + len("pow")
    end = find_matching_paren(text, open_index)
    if end < 0:
      result.append(text[start:])
      break

    args = split_function_args(text[open_index + 1:end])
    if args is None:
      result.append(text[start:end + 1])
      index = end + 1
      continue

    base = strip_outer_parentheses(replace_pow_calls(args[0]))
    exponent = strip_outer_parentheses(replace_pow_calls(args[1]))

    if strip_outer_parentheses(base) == "0.1e1":
      result.append("0.1e1")
    else:
      exponent_key = re.sub(r"\s+", "", exponent)
      replacer = pow_replacements.get(exponent_key)
      if replacer is None:
        result.append("pow({}, {})".format(base, exponent))
      else:
        result.append(replacer(base))

    index = end + 1

  return "".join(result)


def maple2c_replace(text, extra_replace=()):
  '''Performs a series of string replacements in the maple generated C code'''
  
  # The replacements have to be made in order
  math_replace_before_pow = (
    (r"_s_",     r"*"),
    (r"_a_",     r"->"),
    (r"_d_",     r"."),
    (r"_(\d+)_",  r"[\1]"),
    # convert constants like 0.225000000e-1 to 0.225e-1
    (r"0+e", r"e"),
    # convert numerical value of pi to constant
    (r"0.31415926535897932385e1", r"M_PI"),
    # have to do it here, as both Dirac(x) and Dirac(n, x) can appear
    (r"Dirac\(.*?\)", r"0.0"),
    # the derivative of the signum is 0 for us
    (r"signum\(1.*\)", r"0.0"),
  )

  math_replace_after_pow = (
    # cleaning up constant expressions
    (r"sqrt\(0.2e1\)",            r"M_SQRT2"),
    (r"POW_1_3\(0.2e1\)",         r"M_CBRT2"),
    (r"POW_1_3\(0.3e1\)",         r"M_CBRT3"),
    (r"POW_1_3\(0.4e1\)",         r"M_CBRT4"),
    (r"POW_1_3\(0.5e1\)",         r"M_CBRT5"),
    (r"POW_1_3\(0.6e1\)",         r"M_CBRT6"),
    (r"POW_1_3\(M_PI\)",          r"M_CBRTPI"),
  )

  # zk_0_ unfortunatly appears in some expressions
  res = re.search(r"zk_0_ = (.*);", text)
  if res:
    text = re.sub(r"zk_0_(?! =)", "(" + res.group(1) + ")", text)

  # standard replacements
  for str1, str2 in math_replace_before_pow:
    text = re.sub(str1, str2, text)

  text = replace_pow_calls(text)

  for str1, str2 in math_replace_after_pow:
    text = re.sub(str1, str2, text)

  # other specific replacements
  for str1, str2 in extra_replace:
    text = re.sub(str1, str2, text)
  
  return text


def maple_run(params, mtype, code, derivatives, start_order):
  '''Creates the maple file, runs maple, and returns the definition
  of the variables and the c-code
  '''
  
  # open maple file
  from tempfile import mkstemp
  fd, mfilename = mkstemp(suffix=".mpl", text=True)
  fh = os.fdopen(fd, "w")
  
  fh.write('''
Polarization := "{}":
Digits := 20:             (* constants will have 20 digits *)
interface(warnlevel=0):   (* supress all warnings          *)
with(CodeGeneration):

$include "{}"

{}
'''.format(mtype, params["maple_file"], code))
  fh.close()

  # include dirs for maple
  incdirs = ("maple/if_mgga",
             "maple",
             "maple/lda_exc",  "maple/lda_vxc",
             "maple/gga_exc",  "maple/gga_vxc",
             "maple/mgga_exc", "maple/mgga_vxc"
  )
  maple_inc = ["-I" + params["srcdir"] + "/" + i for i in incdirs]

  # run maple
  run = subprocess.run(
    ["maple"] + maple_inc + ["-q", "-u", mfilename], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
  os.remove(mfilename)
  c_code = run.stdout
  if run.returncode != 0 or "Error," in c_code or "Maple initialization failed" in c_code:
    if run.stderr:
      print(run.stderr, file=sys.stderr)
    if c_code:
      print(c_code, file=sys.stderr)
    sys.exit(1)
  if " = " not in c_code:
    print("Maple produced no C assignments", file=sys.stderr)
    sys.exit(1)

  test_1 = ("zk", "vrho", "v2rho2", "v3rho3", "v4rho4", "v5rho5")
  total_order = start_order

  variables  = ["", "", "", "", "", ""]
  n_var = [0, 0, 0, 0, 0, 0]
  new_c_code = ["", "", "", "", "", ""]

  # this adds a new definition of a local variable
  def add_variable(to_add):
    # define 8 variables per line
    if n_var[total_order] % 8 == 0:
      if n_var[total_order] != 0:
        variables[total_order] += ";\n"
      variables[total_order] += "  double "
    else:
      variables[total_order] += ", "
    n_var[total_order] += 1

    variables[total_order] += to_add
  
  # for avoiding compilation when high order derivatives are enabled
  #if start_order != 0:
  #  new_c_code[total_order] += "  if(order < " + str(start_order) + ") return;\n\n\n"

  if params["n_features"] > 0:
    derivative_info = {}
    output_varnames = set()
    last_derivatives = {}
    pending_output_order = None
    pending_output_test = None
    pending_output_lines = []

    for order_index, der_order in enumerate(derivatives):
      actual_order = start_order + order_index
      last_derivatives[actual_order] = der_order[-1][1]
      for der in der_order:
        varname, varorder = parse_output_name(der[1])
        derivative_info[der[1]] = (actual_order, varname, str(varorder))
        output_varnames.add(varname)

    def replace_output_references(line):
      lhs, sep, rhs = line.partition("=")
      if sep == "":
        return line

      rhs = re.sub(
        r"([A-Za-z][A-Za-z0-9]*)_(\d+)_",
        lambda res: "t{}{}".format(res.group(1), res.group(2))
        if res.group(1) in output_varnames else res.group(0),
        rhs)
      return lhs + sep + rhs

    def flush_feature_output_block():
      nonlocal pending_output_order, pending_output_test, pending_output_lines

      if pending_output_test is None:
        return

      new_c_code[pending_output_order] += "  if(" + pending_output_test + ") {\n"
      for output_line in pending_output_lines:
        new_c_code[pending_output_order] += "    " + output_line + "\n"
      new_c_code[pending_output_order] += "  }\n\n"

      pending_output_order = None
      pending_output_test = None
      pending_output_lines = []

    for line in c_code.splitlines():
      line = replace_output_references(line)
      found = False
      res = re.match(r"\s*([A-Za-z][A-Za-z0-9]*_\d+_)\s*=", line)

      if res and res.group(1) in derivative_info:
        out_name = res.group(1)
        actual_order, varname, varorder = derivative_info[out_name]
        assert actual_order == total_order

        add_variable("t" + varname + varorder)
        line = re.sub(r"(\S+)_(\d+)_\s*=\s*(.*);",
                      "t" + varname + varorder + r" = \3;", line)
        new_c_code[total_order] += "  " + line + "\n\n"

        test = "out->" + varname + " != NULL"
        if "lapl" in out_name:
          test += " && (p->info->flags & XC_FLAGS_NEEDS_LAPLACIAN)"
        if "tau" in out_name:
          test += " && (p->info->flags & XC_FLAGS_NEEDS_TAU)"

        test += " && (p->info->flags & XC_FLAGS_HAVE_" + der_name[total_order] + ")"
        output_line = "out->{}[ip*p->dim.{} + {}] += t{}{};".format(
          varname, varname, varorder, varname, varorder)

        if pending_output_test != test or pending_output_order != total_order:
          flush_feature_output_block()
          pending_output_order = total_order
          pending_output_test = test

        pending_output_lines.append(output_line)

        found = True
        if out_name == last_derivatives[total_order]:
          flush_feature_output_block()
          variables[total_order] += ";\n"
          total_order += 1

      if not found:
        res = re.match(r"(t\d+) =", line)
        if res:
          add_variable(res.group(1))

        new_c_code[total_order] += "  " + line + "\n"

    flush_feature_output_block()

  else:
  # we check for strings like 'vrho_0_ = ' and put some
  # relevant if conditions in front
    for line in c_code.splitlines():

      found = False
      # for each order
      for der_order in derivatives:
        # Search the last derivative for each order
        last_derivative = der_order[-1][1]
        # for unpolarized calculation, last derivative is the '0'
        if mtype != "pol" and params["n_features"] == 0:
          last_derivative = re.sub(r"_\d+_", "_0_", last_derivative)

        new_order = re.match(r"\s*" + last_derivative + r"\s*=", line) is not None

        # for each of the derivatives in a given order
        for der in der_order:
          varname  = re.sub(r"_.*", "", der[1])
          varorder = re.sub(r".*_(\d+)_", r"\1", der[1])

          # search for a vrho = statement
          if re.match(r"\s*?" + der[1] + r"\s*=", line):

            if (mtype == "pol") or (params["n_features"] > 0) or (varorder == "0"):
              # we define a new variable (such as tvrho0) to keep the value
              add_variable("t" + varname + varorder)
              line = re.sub(r"(\S+)_(\d+)_\s*=\s*(.*);",
                            "t" + varname + varorder + r" = \3;", line)
              new_c_code[total_order] += "  " + line + "\n\n"

              # build the if clause to assign the variable
              test = "out->" + varname + " != NULL"

              if not re.search(r"lapl", der[1]) is None:
                test += " && (p->info->flags & XC_FLAGS_NEEDS_LAPLACIAN)"

              if not re.search(r"tau", der[1]) is None:
                test += " && (p->info->flags & XC_FLAGS_NEEDS_TAU)"

              test += " && (p->info->flags & XC_FLAGS_HAVE_" + \
                der_name[total_order] + ")"
              new_c_code[total_order] += "  if(" + test + ")\n"

              # add instead of assigning. We are still missing a global constant
              # that can be useful in building hybrid combinations
              new_c_code[total_order] += "    out->{}[ip*p->dim.{} + {}] += t{}{};\n\n".format(varname, varname, varorder, varname, varorder)

            found = True
            break

          # find if vrho_0_ is on the right of the = sign
          # this has necessarily to be defined before the
          # left-hand side of the assignement
          while re.search(r"=.*" + varname + r"_\d+_", line):
            line = re.sub(r"(=.*)" + varname + r"_(\d+)_", r"\1t" + varname + r"\2", line)

        # if last variable of this order increment total_order
        if new_order:
          variables[total_order] += ";\n"
          total_order += 1

      if not found:
        res = re.match(r"(t\d+) =", line)
        if res: add_variable(res.group(1))

        new_c_code[total_order] += "  " + line + "\n"

  # perform the necessary replacements
  for i in range(len(new_c_code)):
    new_c_code[i] = maple2c_replace(new_c_code[i], params["replace"])

  return variables, new_c_code


def maple2c_run(params, variables, derivatives, variants, start_order, input_args, output_args):

  # open file to write to
  fname = params['srcdir'] + "/src/maple2c/" + \
    params['functype']  + "/" + params['functional'] + ".c"
  os.makedirs(os.path.dirname(fname), exist_ok=True)

  from io import StringIO
  out = StringIO()
    
  print_c_header(params, out)

  test_2 = ("EXC", "VXC", "FXC", "KXC", "LXC", "MXC")

  out.write("#define IFXC_MAPLE2C_FLAGS (")
  for i in range(start_order, params['maxorder'] + 1):
    if i != start_order:
      out.write(" | ")
    out.write("IFXC_MGGA_FLAGS_HAVE_" + test_2[i])
  out.write(")\n\n")

  for mtype, variant in variants.items():
    if isinstance(variant, dict) and "batches" in variant:
      batch_specs = variant["batches"]
    elif isinstance(variant, dict):
      batch_specs = [{
        "code": variant["code"],
        "derivatives": variant.get("derivatives", derivatives),
        "start_order": start_order,
      }]
    else:
      batch_specs = [{
        "code": variant,
        "derivatives": derivatives,
        "start_order": start_order,
      }]

    batches = [None] * len(batch_specs)

    def run_batch(batch):
      batch_start = batch.get("start_order", start_order)
      batch_label = batch.get("label")
      if batch_label is not None:
        print("maple2c: running {}".format(batch_label), file=sys.stderr, flush=True)
      vars_def, c_code = maple_run(params, mtype, batch["code"], batch["derivatives"], batch_start)
      if batch_label is not None:
        print("maple2c: finished {}".format(batch_label), file=sys.stderr, flush=True)
      return {
        "start_order": batch_start,
        "end_order": batch_start + len(batch["derivatives"]) - 1,
        "vars_def": vars_def,
        "c_code": c_code,
      }

    if len(batch_specs) > 1:
      requested_workers = params.get("max_workers", os.cpu_count() or 1)
      max_workers = min(len(batch_specs), max(1, requested_workers))
      if max_workers == 1:
        for i, batch in enumerate(batch_specs):
          batches[i] = run_batch(batch)
      else:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
          futures = [executor.submit(run_batch, batch) for batch in batch_specs]
          for i, future in enumerate(futures):
            batches[i] = future.result()
    else:
      batches[0] = run_batch(batch_specs[0])

    for order in range(start_order, params['maxorder'] + 1):
      out.write('''
#ifndef XC_DONT_COMPILE_{}
GPU_DEVICE_FUNCTION static inline void
func_{}_{}(const xc_func_type *p, size_t ip, {}, {})
{{ 
'''.format(der_name[order].upper(),
           der_name[order].lower(), mtype,
           input_args, output_args))

      active_batches = []

      for batch in batches:
        if order < batch["start_order"]:
          continue

        batch_end = min(order, batch["end_order"])
        has_content = params["prefix"] != ""
        for order2 in range(batch["start_order"], batch_end + 1):
          if batch["vars_def"][order2] != "" or batch["c_code"][order2] != "":
            has_content = True

        if not has_content:
          continue

        active_batches.append((batch, batch_end))

      use_scopes = len(active_batches) > 1

      for batch, batch_end in active_batches:
        if use_scopes:
          out.write("  {\n")

        for order2 in range(batch["start_order"], batch_end + 1):
          if batch["vars_def"][order2] != "":
            out.write(batch["vars_def"][order2] + "\n")

        if params["prefix"] != "":
          out.write(params["prefix"] + "\n")

        for order2 in range(batch["start_order"], batch_end + 1):
          out.write(batch["c_code"][order2])

        if use_scopes:
          out.write("  }\n")

      out.write("}\n\n")
      out.write("#endif\n\n")
      
  new_text = out.getvalue()
  out.close()
  with open(fname, "w") as fh:
    fh.write(new_text)
