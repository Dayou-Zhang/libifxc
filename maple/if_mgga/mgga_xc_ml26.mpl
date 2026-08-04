# SPDX-License-Identifier: MPL-2.0

(*
 mgga_xc_ml26.mpl
 IF-MGGA source for the ML26 69-feature semilocal vector.
*)

(* type: if_mgga *)
(* feature_set: ml26 *)
(* nfeatures: 69 *)
(* max_order: 2 *)
(* variables: rho sigma tau *)
(* feature_batches: false *)

(* Reuse the ML25 expression source verbatim, then append the three CS1
   contributions isolated by the development switch parameters. *)
$include "mgga_xc_ml25.mpl"

ifxc_feature_set_key := "ml26":
ifxc_feature_count := 69:

ml26_feature_keys := [
  op(ml25_feature_keys),
  "cs1_same_spin_1",
  "cs1_same_spin_2",
  "cs1_opposite_spin"
]:

ml26_cs1_gamma := 0.006:
ml26_cs1_d := 0.349:
ml26_cs1_C1 := -0.018897:
ml26_cs1_C2 := 0.155240:
ml26_cs1_C3 := -0.159068:
ml26_cs1_C4 := 0.007953:

ml26_cs1_same_spin_prefactor := (rs, z) ->
  opz_pow_n(z, 1)/2 * n_spin(rs, z)^(1/3)
  / (n_spin(rs, z)^(1/3) + ml26_cs1_d):

ml26_cs1_same_spin_constant := (rs, z) ->
  ml26_cs1_same_spin_prefactor(rs, z) * ml26_cs1_C1:

ml26_cs1_same_spin_gradient := (rs, z, xs) ->
  ml26_cs1_same_spin_prefactor(rs, z)
  * ml26_cs1_C2 * ml26_cs1_gamma^2 * xs^4
  / (1 + ml26_cs1_gamma * xs^2)^2:

ml26_cs1_opposite_spin := (rs, z, xt) ->
  (1 - z^2)/4 / (1 + ml26_cs1_d * n_total(rs)^(-1/3))
  * (ml26_cs1_C3
     + ml26_cs1_C4 * ml26_cs1_gamma^2 * xt^4
       / (1 + ml26_cs1_gamma * xt^2)^2):

ml26_feature_exprs := proc(rs, z, xt, xs0, xs1, u0, u1, t0, t1)
  return [
    op(ml25_feature_exprs(rs, z, xt, xs0, xs1, u0, u1, t0, t1)),
    ml26_cs1_same_spin_constant(rs, z)
      + ml26_cs1_same_spin_constant(rs, -z),
    ml26_cs1_same_spin_gradient(rs, z, xs0)
      + ml26_cs1_same_spin_gradient(rs, -z, xs1),
    ml26_cs1_opposite_spin(rs, z, xt)
  ]:
end proc:

f := proc(rs, z, xt, xs0, xs1, u0, u1, t0, t1)
  return ml26_feature_exprs(rs, z, xt, xs0, xs1, u0, u1, t0, t1):
end proc:
