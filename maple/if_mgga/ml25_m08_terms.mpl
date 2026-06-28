(*
 ML25 M08 correlation basis terms used for MN15 A/B features.

 This helper returns the first nine A and B basis terms directly and does not
 assign the generic Maple2C scalar `f`.
*)

$define lda_c_pw_params
$define lda_c_pw_modified_params
$include "lda_c_pw.mpl"

$define gga_c_pbe_params
$include "gga_c_pbe.mpl"

ml25_m08_a_basis := proc(w, f_pw_rs)
  return [
    f_pw_rs,
    w * f_pw_rs,
    w^2 * f_pw_rs,
    w^3 * f_pw_rs,
    w^4 * f_pw_rs,
    w^5 * f_pw_rs,
    w^6 * f_pw_rs,
    w^7 * f_pw_rs,
    w^8 * f_pw_rs
  ]:
end proc:

ml25_m08_b_basis := proc(w, f_pw_rs, f_pbe_rs)
  local delta;

  delta := f_pbe_rs - f_pw_rs:
  return [
    delta,
    w * delta,
    w^2 * delta,
    w^3 * delta,
    w^4 * delta,
    w^5 * delta,
    w^6 * delta,
    w^7 * delta,
    w^8 * delta
  ]:
end proc:

ml25_m08_terms := proc(rs, z, xt, xs0, xs1, t0, t1)
  local t, w, f_pw_rs, f_pbe_rs;

  t := 2^(2/3) * t_total(z, t0, t1):
  w := mgga_w(t):
  f_pw_rs := f_pw(rs, z):
  f_pbe_rs := f_pbe(rs, z, xt, xs0, xs1):

  return [
    op(ml25_m08_a_basis(w, f_pw_rs)),
    op(ml25_m08_b_basis(w, f_pw_rs, f_pbe_rs))
  ]:
end proc:
