(*
 mgga_xc_ml25.mpl
 IF-MGGA source for the ML25 66-feature vector.
*)

(* type: if_mgga *)
(* feature_set: ml25 *)
(* nfeatures: 66 *)
(* max_order: 3 *)
(* variables: rho sigma tau *)
(* feature_batches: false *)

ifxc_feature_set_key := "ml25":
ifxc_feature_count := 66:

ml25_feature_keys := [
  "ml25.lak_x",
  "ml25.lak_c",
  "ml25.lyp_t1",
  "ml25.lyp_t2",
  "ml25.lyp_t3",
  "ml25.lyp_t4",
  "ml25.lyp_t5",
  "ml25.lyp_t6",
  "ml25.mn15_c01",
  "ml25.mn15_c02",
  "ml25.mn15_c03",
  "ml25.mn15_c04",
  "ml25.mn15_c05",
  "ml25.mn15_c06",
  "ml25.mn15_c07",
  "ml25.mn15_c08",
  "ml25.mn15_c09",
  "ml25.mn15_c10",
  "ml25.mn15_c11",
  "ml25.mn15_c12",
  "ml25.mn15_c13",
  "ml25.mn15_c14",
  "ml25.mn15_c15",
  "ml25.mn15_c16",
  "ml25.mn15_c17",
  "ml25.mn15_c18",
  "ml25.mn15_c19",
  "ml25.mn15_c20",
  "ml25.mn15_c21",
  "ml25.mn15_c22",
  "ml25.mn15_c23",
  "ml25.mn15_c24",
  "ml25.mn15_c25",
  "ml25.mn15_c26",
  "ml25.mn15_c27",
  "ml25.mn15_c28",
  "ml25.mn15_c29",
  "ml25.mn15_c30",
  "ml25.mn15_c31",
  "ml25.mn15_c32",
  "ml25.mn15_c33",
  "ml25.mn15_c34",
  "ml25.mn15_c35",
  "ml25.mn15_c36",
  "ml25.mn15_c37",
  "ml25.mn15_c38",
  "ml25.mn15_c39",
  "ml25.mn15_c40",
  "ml25.mn15_a01",
  "ml25.mn15_a02",
  "ml25.mn15_a03",
  "ml25.mn15_a04",
  "ml25.mn15_a05",
  "ml25.mn15_a06",
  "ml25.mn15_a07",
  "ml25.mn15_a08",
  "ml25.mn15_a09",
  "ml25.mn15_b01",
  "ml25.mn15_b02",
  "ml25.mn15_b03",
  "ml25.mn15_b04",
  "ml25.mn15_b05",
  "ml25.mn15_b06",
  "ml25.mn15_b07",
  "ml25.mn15_b08",
  "ml25.mn15_b09"
]:

$include "ml25_lak.mpl"
$include "ml25_lyp.mpl"
$include "ml25_mn12_terms.mpl"
$include "ml25_m08_terms.mpl"

ml25_feature_exprs := proc(rs, z, xt, xs0, xs1, u0, u1, t0, t1)
  local lyp_terms, mn12_terms, m08_terms;

  lyp_terms := ml25_lyp_terms(rs, z, xt, xs0, xs1):
  mn12_terms := ml25_mn12_c_terms(rs, z, xs0, xs1, u0, u1, t0, t1):
  m08_terms := ml25_m08_terms(rs, z, xt, xs0, xs1, t0, t1):

  return [
    ml25_lak_x_f(rs, z, xt, xs0, xs1, u0, u1, t0, t1),
    ml25_lak_c_f(rs, z, xt, xs0, xs1, u0, u1, t0, t1),
    op(lyp_terms),
    op(mn12_terms),
    op(m08_terms)
  ]:
end proc:

f := proc(rs, z, xt, xs0, xs1, u0, u1, t0, t1)
  return ml25_feature_exprs(rs, z, xt, xs0, xs1, u0, u1, t0, t1):
end proc:
