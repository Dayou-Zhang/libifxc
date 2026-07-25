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
  "ml25.mn15_cc000",
  "ml25.mn15_cc001",
  "ml25.mn15_cc002",
  "ml25.mn15_cc003",
  "ml25.mn15_cc004",
  "ml25.mn15_cc005",
  "ml25.mn15_cc010",
  "ml25.mn15_cc011",
  "ml25.mn15_cc012",
  "ml25.mn15_cc013",
  "ml25.mn15_cc014",
  "ml25.mn15_cc020",
  "ml25.mn15_cc021",
  "ml25.mn15_cc022",
  "ml25.mn15_cc023",
  "ml25.mn15_cc030",
  "ml25.mn15_cc031",
  "ml25.mn15_cc032",
  "ml25.mn15_cc100",
  "ml25.mn15_cc101",
  "ml25.mn15_cc102",
  "ml25.mn15_cc103",
  "ml25.mn15_cc104",
  "ml25.mn15_cc110",
  "ml25.mn15_cc111",
  "ml25.mn15_cc112",
  "ml25.mn15_cc113",
  "ml25.mn15_cc120",
  "ml25.mn15_cc121",
  "ml25.mn15_cc122",
  "ml25.mn15_cc200",
  "ml25.mn15_cc201",
  "ml25.mn15_cc202",
  "ml25.mn15_cc203",
  "ml25.mn15_cc210",
  "ml25.mn15_cc211",
  "ml25.mn15_cc212",
  "ml25.mn15_cc300",
  "ml25.mn15_cc301",
  "ml25.mn15_cc302",
  "ml25.mn15_a0",
  "ml25.mn15_a1",
  "ml25.mn15_a2",
  "ml25.mn15_a3",
  "ml25.mn15_a4",
  "ml25.mn15_a5",
  "ml25.mn15_a6",
  "ml25.mn15_a7",
  "ml25.mn15_a8",
  "ml25.mn15_b0",
  "ml25.mn15_b1",
  "ml25.mn15_b2",
  "ml25.mn15_b3",
  "ml25.mn15_b4",
  "ml25.mn15_b5",
  "ml25.mn15_b6",
  "ml25.mn15_b7",
  "ml25.mn15_b8"
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
