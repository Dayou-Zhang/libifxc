# SPDX-License-Identifier: MPL-2.0

(*
 ML25 LYP correlation feature terms.

 This helper defines the six ML25 LYP basis expressions in one list and does
 not assign the generic Maple2C scalar `f`.
*)

$include "util.mpl"

ml25_lyp_params_a := 0.04918:
ml25_lyp_params_b := 0.132:
ml25_lyp_params_c := 0.2533:
ml25_lyp_params_d := 0.349:

ml25_lyp_Cf := 3/10 * (3*Pi^2)^(2/3):

ml25_lyp_omega := rr ->
  ml25_lyp_params_b*exp(-ml25_lyp_params_c*rr)/(1 + ml25_lyp_params_d*rr):

ml25_lyp_delta := rr ->
  (ml25_lyp_params_c + ml25_lyp_params_d/(1 + ml25_lyp_params_d*rr))*rr:

ml25_lyp_aux6 := 1/2^(8/3):
ml25_lyp_aux4 := ml25_lyp_aux6/4:
ml25_lyp_aux5 := ml25_lyp_aux4/(9*2):

ml25_lyp_t1 := (rr, z) ->
  -(1 - z^2)/(1 + ml25_lyp_params_d*rr):

ml25_lyp_t2 := (rr, z, xt) ->
  -xt^2*((1 - z^2)*(47 - 7*ml25_lyp_delta(rr))/(4*18) - 2/3):

ml25_lyp_t3 := z ->
  -ml25_lyp_Cf/2*(1 - z^2)*(opz_pow_n(z, 8/3) + opz_pow_n(-z, 8/3)):

ml25_lyp_t4 := (rr, z, xs0, xs1) ->
  ml25_lyp_aux4*(1 - z^2)*(5/2 - ml25_lyp_delta(rr)/18)*
  (xs0^2*opz_pow_n(z, 8/3) + xs1^2*opz_pow_n(-z, 8/3)):

ml25_lyp_t5 := (rr, z, xs0, xs1) ->
  ml25_lyp_aux5*(1 - z^2)*(ml25_lyp_delta(rr) - 11)*
  (xs0^2*opz_pow_n(z, 11/3) + xs1^2*opz_pow_n(-z, 11/3)):

ml25_lyp_t6 := (z, xs0, xs1) ->
  -ml25_lyp_aux6*(2/3*(xs0^2*opz_pow_n(z, 8/3) + xs1^2*opz_pow_n(-z, 8/3))
  - opz_pow_n(z, 2)*xs1^2*opz_pow_n(-z, 8/3)/4
  - opz_pow_n(-z, 2)*xs0^2*opz_pow_n(z, 8/3)/4):

ml25_lyp_terms_rr := proc(rr, z, xt, xs0, xs1)
  return [
    ml25_lyp_params_a * ml25_lyp_t1(rr, z),
    ml25_lyp_params_a * ml25_lyp_omega(rr) * ml25_lyp_t2(rr, z, xt),
    ml25_lyp_params_a * ml25_lyp_omega(rr) * ml25_lyp_t3(z),
    ml25_lyp_params_a * ml25_lyp_omega(rr) * ml25_lyp_t4(rr, z, xs0, xs1),
    ml25_lyp_params_a * ml25_lyp_omega(rr) * ml25_lyp_t5(rr, z, xs0, xs1),
    ml25_lyp_params_a * ml25_lyp_omega(rr) * ml25_lyp_t6(z, xs0, xs1)
  ]:
end proc:

ml25_lyp_terms := (rs, z, xt, xs0, xs1) ->
  ml25_lyp_terms_rr(rs/RS_FACTOR, z, xt, xs0, xs1):
