(*
 ML25 MN12/MN15 C exchange basis terms.

 This helper returns the 40 basis terms in the same order as the old
 one-coefficient wrappers and does not assign the generic Maple2C scalar `f`.
*)

$include "util.mpl"

ml25_mn12_omega_x := 2.5:
ml25_mn12_gamma_x := 0.004:

ml25_mn12_vx := (rs, z) ->
  1/(1 + rs/(ml25_mn12_omega_x*RS_FACTOR)*(2/(1 + z))^(1/3)):

ml25_mn12_ux := x ->
  ml25_mn12_gamma_x*x^2/(1 + ml25_mn12_gamma_x*x^2):

ml25_mn12_wx := t ->
  (K_FACTOR_C - t)/(K_FACTOR_C + t):

ml25_mn12_c_basis := proc(rs, z, x, u, t)
  local ux, ux2, ux3, vx, vx2, vx3, w;

  ux := ml25_mn12_ux(x):
  ux2 := ux*ux:
  ux3 := ux2*ux:
  vx := ml25_mn12_vx(rs, z):
  vx2 := vx*vx:
  vx3 := vx2*vx:
  w := ml25_mn12_wx(t):

  return [
    1, w, w^2, w^3, w^4, w^5,
    ux, ux*w, ux*w^2, ux*w^3, ux*w^4,
    ux2, ux2*w, ux2*w^2, ux2*w^3,
    ux3, ux3*w, ux3*w^2,
    vx, vx*w, vx*w^2, vx*w^3, vx*w^4,
    ux*vx, ux*vx*w, ux*vx*w^2, ux*vx*w^3,
    ux2*vx, ux2*vx*w, ux2*vx*w^2,
    vx2, vx2*w, vx2*w^2, vx2*w^3,
    ux*vx2, ux*vx2*w, ux*vx2*w^2,
    vx3, vx3*w, vx3*w^2
  ]:
end proc:

ml25_mn12_c_terms := proc(rs, z, xs0, xs1, u0, u1, t0, t1)
  local i, up_terms, down_terms;

  up_terms := ml25_mn12_c_basis(rs, z_thr(z), xs0, u0, t0):
  down_terms := ml25_mn12_c_basis(rs, z_thr(-z), xs1, u1, t1):

  return [
    seq(
      my_piecewise3(screen_dens(rs, z), 0, lda_x_spin(rs, z_thr(z))*op(i, up_terms)) +
      my_piecewise3(screen_dens(rs, -z), 0, lda_x_spin(rs, z_thr(-z))*op(i, down_terms)),
      i = 1..40)
  ]:
end proc:
