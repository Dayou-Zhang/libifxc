(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

$include "gga_c_lyp.mpl"

params_a_a := 0.04918:
params_a_b := 0.132:
params_a_c := 0.2533:
params_a_d := 0.349:

f_lyp_rr := (rr, z, xt, xs0, xs1) -> params_a_a*lyp_omega(rr)*(
  lyp_t4(rr, z, xs0, xs1)
):

(* rr = rs/RS_FACTOR is equal to n_total(rs)^(-1/3) *)
f_lyp := (rs, z, xt, xs0, xs1) -> f_lyp_rr(rs/RS_FACTOR, z, xt, xs0, xs1):

f  := (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) -> f_lyp(rs, z, xt, xs0, xs1):
