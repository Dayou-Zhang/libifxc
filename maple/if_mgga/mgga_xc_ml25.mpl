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
  "lak_x",
  "lak_c",
  "lyp_t1",
  "lyp_t2",
  "lyp_t3",
  "lyp_t4",
  "lyp_t5",
  "lyp_t6",
  "mn15_cc000",
  "mn15_cc001",
  "mn15_cc002",
  "mn15_cc003",
  "mn15_cc004",
  "mn15_cc005",
  "mn15_cc010",
  "mn15_cc011",
  "mn15_cc012",
  "mn15_cc013",
  "mn15_cc014",
  "mn15_cc020",
  "mn15_cc021",
  "mn15_cc022",
  "mn15_cc023",
  "mn15_cc030",
  "mn15_cc031",
  "mn15_cc032",
  "mn15_cc100",
  "mn15_cc101",
  "mn15_cc102",
  "mn15_cc103",
  "mn15_cc104",
  "mn15_cc110",
  "mn15_cc111",
  "mn15_cc112",
  "mn15_cc113",
  "mn15_cc120",
  "mn15_cc121",
  "mn15_cc122",
  "mn15_cc200",
  "mn15_cc201",
  "mn15_cc202",
  "mn15_cc203",
  "mn15_cc210",
  "mn15_cc211",
  "mn15_cc212",
  "mn15_cc300",
  "mn15_cc301",
  "mn15_cc302",
  "mn15_a0",
  "mn15_a1",
  "mn15_a2",
  "mn15_a3",
  "mn15_a4",
  "mn15_a5",
  "mn15_a6",
  "mn15_a7",
  "mn15_a8",
  "mn15_b0",
  "mn15_b1",
  "mn15_b2",
  "mn15_b3",
  "mn15_b4",
  "mn15_b5",
  "mn15_b6",
  "mn15_b7",
  "mn15_b8"
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
