(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_c_m08_params *params;

  assert(p->params != NULL);
  params = (mgga_c_m08_params * )(p->params);
*)


$define lda_c_pw_params
$define lda_c_pw_modified_params
$include "lda_c_pw.mpl"

$define gga_c_pbe_params
$include "gga_c_pbe.mpl"

(* the prefactor of t was chosen to get the right K_FACTOR_C in mgga_series_from_w *)
m08_f := proc(rs, z, xt, xs0, xs1, ts0, ts1)
  local t, w, f_pw_rs, f_pbe_rs;
  t := 2^(2/3) * t_total(z, ts0, ts1);
  w := mgga_w(t);
  f_pw_rs := f_pw(rs, z);
  f_pbe_rs := f_pbe(rs, z, xt, xs0, xs1);
  return
    + mgga_series_from_w(params_a_m08_a, 12, w) * f_pw_rs
    + mgga_series_from_w(params_a_m08_b, 12, w) * (f_pbe_rs - f_pw_rs);
end proc:

f := (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) ->
  m08_f(rs, z, xt, xs0, xs1, ts0, ts1):
