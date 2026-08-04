# SPDX-License-Identifier: MPL-2.0

(*
 ML25 LAK exchange and correlation feature terms.

 This helper keeps LAK exchange and correlation in one source file and does
 not assign the generic Maple2C scalar `f`.
*)

$include "util.mpl"

(* Parameters from page 3. *)
lak_h0x := 1.174:
lak_mu_ax := -(97 + 3*lak_h0x + sqrt((3*lak_h0x)^2 + 74166*lak_h0x - 64175))/1200:
lak_nu_a := (73-50*lak_mu_ax)/5000:
lak_mu_sx := (10+60*lak_mu_ax)/81:
lak_nu_s := -(1606 - 50*lak_mu_ax)/18225:
lak_bx := 4.9479:
lak_ax := 1.1:
lak_anum := 5:

(* SI eq 2. *)
lak_fsa := (s, a) ->
  lak_h0x*lak_gx(s) + (1-lak_fx(a))*(lak_h1x(s)-lak_h0x)*lak_gnum(s):

(* SI eq 4. *)
lak_gx := s -> 1 - m_recexp(sqrt(s)/lak_bx):

(* SI eq 5. *)
lak_fx0 := a -> 2/Pi * arctan(Pi/2*(lak_c1*(a-1)/a + lak_c2*(a-1)^2)):
lak_fx_taylor := a -> eval(convert(eval(series(lak_fx0(b), b=0, 9), csgn=-1), polynom), b=a):
lak_fx := a -> my_piecewise3(a <= DBL_EPSILON, lak_fx_taylor(a), lak_fx0(m_max(a, DBL_EPSILON))):

(* SI eq 6. *)
lak_c1 := lak_mu_ax/(lak_h0x-1):

(* SI eq 7. *)
lak_c2 := (lak_mu_ax + lak_nu_a)/(lak_h0x-1):

(* SI eq 8. *)
lak_h1x := s -> lak_hx_ge4(s) + lak_kx(s)*(lak_ax - lak_hx_ge4(s)):

(* SI eq 9. *)
lak_hx_ge4 := s -> 1 + lak_mu_sx*s^2 + lak_nu_s*s^4 + lak_h0x*(1-lak_gx(s)):

(* SI eq 10. *)
lak_kx := s -> m_recexp((s/lak_ax)^2 * (1+s^2)):

(* SI eq 11. *)
lak_gnum := s -> 1 - m_recexp((s/lak_anum)^2):

(* Exchange feature. *)
ml25_lak_alpha_x := (x, t) -> (t/K_FACTOR_C) * m_max(1 - x^2/(8*t), 1e-10):
ml25_lak_f_x := (x, u, t) -> lak_fsa(x*X2S, ml25_lak_alpha_x(x, t)):

ml25_lak_x_f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(ml25_lak_f_x, rs, z, xs0, xs1, u0, u1, t0, t1):

(* Correlation feature. *)
$define lda_c_pw_params
$define lda_c_pw_modified_params
$include "lda_c_pw.mpl"

(* Parameters from page 4 of SI. *)
lak_mu_a := -lak_mu_ax/2:
lak_b1c := 0.0468:
lak_b2c := 0.205601:
lak_chi_zero := 1.55344:
lak_b3c := 2.85:
lak_ac := 10:

(* SI eq 3. *)
lak_Ax := -3/4 * (3/Pi)^(1/3):

(* SI eq 12. *)
lak_ec_rszta := (rs, z, xt, a) ->
  lak_e0(rs, z, XT2S*xt) +
  (1-lak_fc(rs, a)) * (lak_e1(rs, z, tt(rs, z, xt)) - lak_e0(rs, z, XT2S*xt)) *
  lak_gnum(XT2S*xt):

(* SI eq 13. *)
lak_e0 := (rs, z, s) -> (lak_eclda0(rs) + lak_H0(rs, s))*lak_Gc(z):

(* SI eq 14. *)
lak_eclda0 := rs -> -lak_b1c / (1 + lak_b2c*rs):

(* SI eq 15. *)
lak_H0 := (rs, s) ->
  lak_b1c*log(m_max(DBL_MIN, 1 + lak_w0(rs)*(1 - (1 + 4*lak_chi_zero*s^2)^(-1/4)))):

(* SI eq 16. *)
lak_w0 := rs -> exp(-lak_eclda0(rs)/lak_b1c) - 1:

(* SI eq 17. *)
lak_G_cnst := 2.3631:
lak_Gc := z -> (1 - lak_G_cnst*(2^(1/3) - 1)*f_zeta(z))*(1 - z^12):

(* SI eq 19. *)
lak_fc0 := x -> 2/Pi * arctan(Pi/2 * x):
lak_fc_taylor := x -> eval(convert(eval(series(lak_fc0(y), y=infinity, 9), csgn=1), polynom), y=x):
lak_fc_piecewise := x ->
  my_piecewise3(x >= 1/DBL_EPSILON,
    lak_fc_taylor(m_max(x, 1/DBL_EPSILON)),
    lak_fc0(m_min(x, 1/DBL_EPSILON))):
lak_fc := (rs, a) ->
  my_piecewise3(a < 0,
    2+lak_fc_piecewise(lak_fcGE2(rs)*m_abs(lak_at(rs, a))),
    lak_fc_piecewise(lak_fcGE2(rs)*lak_at(rs, a))):

(* SI eq 20. *)
lak_at := (rs, a) -> (a-1)/(m_max(rs*a, DBL_EPSILON)):

(* SI eq 21. *)
lak_fcGE2 := rs -> lak_beta_at(rs) / (m_min(f_pw(rs, 0)-lak_eclda0(rs), -DBL_EPSILON)):

(* SI eq 22. *)
lak_beta_at := rs -> -lak_Ac * lak_mu_ac(rs):

(* SI eq 23. *)
lak_mu_ac := rs -> lak_Cma(rs)*lak_mu_a - lak_mu_ax:

(* SI eq 24. *)
lak_e1 := (rs, z, t) -> f_pw(rs, z) + lak_H1(rs, z, t):

(* SI eq 25. *)
lak_H1 := (rs, z, t) -> my_piecewise3(lak_beta_t(rs) >= 0, lak_H1p(rs, z, t), lak_H1m(rs, z, t)):

(* SI eq 27. *)
lak_ct := z -> (3*Pi^2/16)^(2/3) / mphi(z)^2:

(* SI eq 29. *)
lak_beta_t := rs -> lak_Ac / lak_ct(0) * lak_mu_sc(rs):

(* SI eq 30. *)
lak_Ac := RS_FACTOR * lak_Ax:

(* SI eq 31. *)
lak_mu_sc := rs -> MU_GE * (lak_Cma(rs)*(1+6*lak_mu_a) - (1+6*lak_mu_ax)):

(* SI eq 33. *)
lak_H1p := (rs, z, t) ->
  lak_gamma*mphi(z)^3*log(m_max(1 + lak_w1(rs, z) *
  (1 - lak_g1p(rs, z, t)) * (1 - lak_g2(rs, z, t) + lak_g3p(rs, z, t)), DBL_MIN)):

(* SI eq 34. *)
lak_H1m := (rs, z, t) ->
  lak_gamma*mphi(z)^3*log(m_max(1 + lak_w1(rs, z) *
  (1 - lak_g1m(rs, z, t)) * (1 - lak_g2(rs, z, t) - lak_g3m(rs, z, t)), DBL_MIN)):

(* SI eq 35. *)
lak_gamma := (1 - log(2))/Pi^2:

(* SI eq 36. *)
lak_w1 := (rs, z) -> exp(-f_pw(rs, z)/(lak_gamma*mphi(z)^3)) - 1:

(* SI eq 37. *)
lak_g1p := (rs, z, t) -> (m_max(1 + 4*lak_A(rs, z)*t^2, DBL_EPSILON))^(-1/4):
lak_g1m := (rs, z, t) -> (m_max(1 - 4*lak_A(rs, z)*t^2, DBL_EPSILON))^(-1/4):

(* SI eq 38. *)
lak_g2 := (rs, z, t) -> (1 + (lak_A(rs, z)*t^2)^2)^(-1):

(* SI eq 39. *)
lak_g3p := (rs, z, t) -> (m_max(1 + lak_ac*lak_A(rs, z)*t^2, DBL_EPSILON))^(-1):

(* SI eq 40. *)
lak_g3m := (rs, z, t) ->
  (m_max(1 - (lak_w1(rs, z)+lak_b3c) * lak_A(rs, z) * t^2, DBL_EPSILON))^(-1):

(* SI eq 41. *)
lak_A := (rs, z) -> lak_beta_t(rs) / (lak_gamma*lak_w1(rs, z)):

(* SI eq 62. *)
lak_Cma := rs ->
  lak_Cs0 * (1 + 0.1*rs^0.65) / ((1 + 0.065*rs^0.9)*(1 + 0.03*rs^1.2)):

(* Constant in SI main text after eq 62. *)
lak_Cs0 := -16*Pi*(3*Pi^2)^(1/3)*2.568/(3000*MU_GE):

ml25_lak_alpha_c := (z, xt, ts0, ts1) ->
  (t_total(z, ts0, ts1) - xt^2/8)/(K_FACTOR_C * t_total(z, 1, 1)):

ml25_lak_ec := (rs, z, xt, xs0, xs1, ts0, ts1) ->
  lak_ec_rszta(rs, z, xt, ml25_lak_alpha_c(z, xt, ts0, ts1)):

ml25_lak_c_f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  ml25_lak_ec(rs, z, xt, xs0, xs1, t0, t1):
