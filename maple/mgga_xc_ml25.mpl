(*
 mgga_xc_ml25.mpl
 Source of truth for libifxc ML25 integral-feature kernels.
 Feature order must match src/features/ifxc_ml25_features.def.
*)

ifxc_feature_set_key := "ml25":
ifxc_feature_count := 66:

$include "mgga_exc/mgga_xc_ml25_common.mpl"

# 0: LAK exchange
$include "mgga_exc/mgga_x_lak.mpl"
lak_x_h := f:

# 1: LAK correlation
$include "mgga_exc/mgga_c_lak.mpl"
lak_c_h := f:

# 2-7: LYP correlation terms
$include "mgga_exc/mgga_xc_ml25_lyp_t1.mpl"
lyp_t1_h := f:
$include "mgga_exc/mgga_xc_ml25_lyp_t2.mpl"
lyp_t2_h := f:
$include "mgga_exc/mgga_xc_ml25_lyp_t3.mpl"
lyp_t3_h := f:
$include "mgga_exc/mgga_xc_ml25_lyp_t4.mpl"
lyp_t4_h := f:
$include "mgga_exc/mgga_xc_ml25_lyp_t5.mpl"
lyp_t5_h := f:
$include "mgga_exc/mgga_xc_ml25_lyp_t6.mpl"
lyp_t6_h := f:

# 8-47: MN15 c terms
$include "mgga_exc/mgga_xc_ml25_mn15_c01.mpl"
mn15_c01_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c02.mpl"
mn15_c02_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c03.mpl"
mn15_c03_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c04.mpl"
mn15_c04_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c05.mpl"
mn15_c05_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c06.mpl"
mn15_c06_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c07.mpl"
mn15_c07_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c08.mpl"
mn15_c08_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c09.mpl"
mn15_c09_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c10.mpl"
mn15_c10_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c11.mpl"
mn15_c11_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c12.mpl"
mn15_c12_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c13.mpl"
mn15_c13_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c14.mpl"
mn15_c14_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c15.mpl"
mn15_c15_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c16.mpl"
mn15_c16_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c17.mpl"
mn15_c17_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c18.mpl"
mn15_c18_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c19.mpl"
mn15_c19_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c20.mpl"
mn15_c20_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c21.mpl"
mn15_c21_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c22.mpl"
mn15_c22_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c23.mpl"
mn15_c23_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c24.mpl"
mn15_c24_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c25.mpl"
mn15_c25_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c26.mpl"
mn15_c26_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c27.mpl"
mn15_c27_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c28.mpl"
mn15_c28_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c29.mpl"
mn15_c29_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c30.mpl"
mn15_c30_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c31.mpl"
mn15_c31_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c32.mpl"
mn15_c32_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c33.mpl"
mn15_c33_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c34.mpl"
mn15_c34_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c35.mpl"
mn15_c35_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c36.mpl"
mn15_c36_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c37.mpl"
mn15_c37_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c38.mpl"
mn15_c38_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c39.mpl"
mn15_c39_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_c40.mpl"
mn15_c40_h := f:

# 48-56: MN15 a terms
$include "mgga_exc/mgga_xc_ml25_mn15_a01.mpl"
mn15_a01_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_a02.mpl"
mn15_a02_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_a03.mpl"
mn15_a03_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_a04.mpl"
mn15_a04_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_a05.mpl"
mn15_a05_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_a06.mpl"
mn15_a06_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_a07.mpl"
mn15_a07_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_a08.mpl"
mn15_a08_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_a09.mpl"
mn15_a09_h := f:

# 57-65: MN15 b terms
$include "mgga_exc/mgga_xc_ml25_mn15_b01.mpl"
mn15_b01_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_b02.mpl"
mn15_b02_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_b03.mpl"
mn15_b03_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_b04.mpl"
mn15_b04_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_b05.mpl"
mn15_b05_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_b06.mpl"
mn15_b06_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_b07.mpl"
mn15_b07_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_b08.mpl"
mn15_b08_h := f:
$include "mgga_exc/mgga_xc_ml25_mn15_b09.mpl"
mn15_b09_h := f:

ifxc_ml25_features := [
  lak_x_h,    # 0, ml25.lak_x
  lak_c_h,    # 1, ml25.lak_c
  lyp_t1_h,   # 2, ml25.lyp_t1
  lyp_t2_h,   # 3, ml25.lyp_t2
  lyp_t3_h,   # 4, ml25.lyp_t3
  lyp_t4_h,   # 5, ml25.lyp_t4
  lyp_t5_h,   # 6, ml25.lyp_t5
  lyp_t6_h,   # 7, ml25.lyp_t6
  mn15_c01_h, # 8, ml25.mn15_c01
  mn15_c02_h, # 9, ml25.mn15_c02
  mn15_c03_h, # 10, ml25.mn15_c03
  mn15_c04_h, # 11, ml25.mn15_c04
  mn15_c05_h, # 12, ml25.mn15_c05
  mn15_c06_h, # 13, ml25.mn15_c06
  mn15_c07_h, # 14, ml25.mn15_c07
  mn15_c08_h, # 15, ml25.mn15_c08
  mn15_c09_h, # 16, ml25.mn15_c09
  mn15_c10_h, # 17, ml25.mn15_c10
  mn15_c11_h, # 18, ml25.mn15_c11
  mn15_c12_h, # 19, ml25.mn15_c12
  mn15_c13_h, # 20, ml25.mn15_c13
  mn15_c14_h, # 21, ml25.mn15_c14
  mn15_c15_h, # 22, ml25.mn15_c15
  mn15_c16_h, # 23, ml25.mn15_c16
  mn15_c17_h, # 24, ml25.mn15_c17
  mn15_c18_h, # 25, ml25.mn15_c18
  mn15_c19_h, # 26, ml25.mn15_c19
  mn15_c20_h, # 27, ml25.mn15_c20
  mn15_c21_h, # 28, ml25.mn15_c21
  mn15_c22_h, # 29, ml25.mn15_c22
  mn15_c23_h, # 30, ml25.mn15_c23
  mn15_c24_h, # 31, ml25.mn15_c24
  mn15_c25_h, # 32, ml25.mn15_c25
  mn15_c26_h, # 33, ml25.mn15_c26
  mn15_c27_h, # 34, ml25.mn15_c27
  mn15_c28_h, # 35, ml25.mn15_c28
  mn15_c29_h, # 36, ml25.mn15_c29
  mn15_c30_h, # 37, ml25.mn15_c30
  mn15_c31_h, # 38, ml25.mn15_c31
  mn15_c32_h, # 39, ml25.mn15_c32
  mn15_c33_h, # 40, ml25.mn15_c33
  mn15_c34_h, # 41, ml25.mn15_c34
  mn15_c35_h, # 42, ml25.mn15_c35
  mn15_c36_h, # 43, ml25.mn15_c36
  mn15_c37_h, # 44, ml25.mn15_c37
  mn15_c38_h, # 45, ml25.mn15_c38
  mn15_c39_h, # 46, ml25.mn15_c39
  mn15_c40_h, # 47, ml25.mn15_c40
  mn15_a01_h, # 48, ml25.mn15_a01
  mn15_a02_h, # 49, ml25.mn15_a02
  mn15_a03_h, # 50, ml25.mn15_a03
  mn15_a04_h, # 51, ml25.mn15_a04
  mn15_a05_h, # 52, ml25.mn15_a05
  mn15_a06_h, # 53, ml25.mn15_a06
  mn15_a07_h, # 54, ml25.mn15_a07
  mn15_a08_h, # 55, ml25.mn15_a08
  mn15_a09_h, # 56, ml25.mn15_a09
  mn15_b01_h, # 57, ml25.mn15_b01
  mn15_b02_h, # 58, ml25.mn15_b02
  mn15_b03_h, # 59, ml25.mn15_b03
  mn15_b04_h, # 60, ml25.mn15_b04
  mn15_b05_h, # 61, ml25.mn15_b05
  mn15_b06_h, # 62, ml25.mn15_b06
  mn15_b07_h, # 63, ml25.mn15_b07
  mn15_b08_h, # 64, ml25.mn15_b08
  mn15_b09_h  # 65, ml25.mn15_b09
]:
