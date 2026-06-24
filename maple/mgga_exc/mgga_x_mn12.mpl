(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_x_mn12_params *params;

  assert(p->params != NULL);
  params = (mgga_x_mn12_params * ) (p->params);
*)

$define lda_x_params
$include "lda_x.mpl"

mn12_omega_x := 2.5:
mn12_gamma_x := 0.004:

mn12_vx := (rs, z) -> 1/(1 + rs/(mn12_omega_x*RS_FACTOR)*(2/(1 + z))^(1/3)):
mn12_ux := x -> mn12_gamma_x*x^2/(1 + mn12_gamma_x*x^2):
mn12_wx := t -> (K_FACTOR_C - t)/(K_FACTOR_C + t):

mn12_poly := proc(first, last, w)
  local i, acc;
  acc := params_a_c[last];
  for i from last - 1 by -1 to first do
    acc := params_a_c[i] + w*acc;
  end do;
  return acc;
end proc:

mn12_pol1 := proc(w)
  return mn12_poly(1, 6, w);
end proc:
mn12_pol2 := proc(w)
  return mn12_poly(7, 11, w);
end proc:
mn12_pol3 := proc(w)
  return mn12_poly(12, 15, w);
end proc:
mn12_pol4 := proc(w)
  return mn12_poly(16, 18, w);
end proc:
mn12_pol5 := proc(w)
  return mn12_poly(19, 23, w);
end proc:
mn12_pol6 := proc(w)
  return mn12_poly(24, 27, w);
end proc:
mn12_pol7 := proc(w)
  return mn12_poly(28, 30, w);
end proc:
mn12_pol8 := proc(w)
  return mn12_poly(31, 34, w);
end proc:
mn12_pol9 := proc(w)
  return mn12_poly(35, 37, w);
end proc:
mn12_pol10 := proc(w)
  return mn12_poly(38, 40, w);
end proc:

mn12_f := proc(rs, z, x, u, t)
  local vx, vx2, vx3, ux, ux2, ux3, w;
  vx := mn12_vx(rs, z);
  vx2 := vx*vx;
  vx3 := vx2*vx;
  ux := mn12_ux(x);
  ux2 := ux*ux;
  ux3 := ux2*ux;
  w := mn12_wx(t);
  return
    + mn12_pol1(w)
    + mn12_pol2(w)*ux
    + mn12_pol3(w)*ux2
    + mn12_pol4(w)*ux3
    + mn12_pol5(w)*vx
    + mn12_pol6(w)*ux*vx
    + mn12_pol7(w)*ux2*vx
    + mn12_pol8(w)*vx2
    + mn12_pol9(w)*ux*vx2
    + mn12_pol10(w)*vx3;
end proc:

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange_nsp(mn12_f, rs, z, xs0, xs1, u0, u1, t0, t1):
