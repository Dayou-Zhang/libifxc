(*
 Copyright (C) 2024 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

(* Build the functional *)
lak_alpha := (x, t) -> (t/K_FACTOR_C) * m_max(1 - x^2/(8*t), 1e-10):
lak_f   := (x, u, t) -> lak_fsa(x*X2S, lak_alpha(x,t)):
f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) -> mgga_exchange(lak_f, rs, z, xs0, xs1, u0, u1, t0, t1):
