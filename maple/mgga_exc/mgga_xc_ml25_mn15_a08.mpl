(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

params_a_m08_a := Array(1..12, 0):
params_a_m08_b := Array(1..12, 0):
params_a_m08_a[8] := 1:

$include "mgga_c_m08.mpl"
