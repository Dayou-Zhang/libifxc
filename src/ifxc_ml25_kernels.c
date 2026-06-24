/*
 Copyright (C) 2016 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

/* ID 0 gives the original MN15 */
#define XC_PREFIX           mn15_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15.c"
#include "mgga_xc_ml25_inc.c"


/* Now include work_mgga functions of each integral */
/* Each work_mgga function is prefixed by macro */

/* Integral 1: LAK exchange */
#define XC_PREFIX           lak_x_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_x_lak.c"
#include "mgga_xc_ml25_inc.c"

/* Integral 2: LAK correlation */
#define XC_PREFIX           lak_c_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_c_lak.c"
#include "mgga_xc_ml25_inc.c"

/* Integral 3-8: LYP correlation (6 terms) */
#define XC_PREFIX           lyp_c_t1_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_lyp_t1.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           lyp_c_t2_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_lyp_t2.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           lyp_c_t3_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_lyp_t3.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           lyp_c_t4_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_lyp_t4.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           lyp_c_t5_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_lyp_t5.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           lyp_c_t6_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_lyp_t6.c"
#include "mgga_xc_ml25_inc.c"

/* Integral 9-48: MN15 c terms (40 terms) */
#define XC_PREFIX           mn15_c01_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c01.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c02_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c02.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c03_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c03.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c04_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c04.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c05_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c05.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c06_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c06.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c07_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c07.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c08_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c08.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c09_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c09.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c10_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c10.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c11_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c11.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c12_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c12.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c13_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c13.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c14_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c14.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c15_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c15.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c16_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c16.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c17_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c17.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c18_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c18.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c19_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c19.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c20_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c20.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c21_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c21.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c22_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c22.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c23_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c23.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c24_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c24.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c25_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c25.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c26_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c26.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c27_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c27.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c28_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c28.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c29_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c29.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c30_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c30.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c31_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c31.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c32_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c32.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c33_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c33.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c34_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c34.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c35_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c35.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c36_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c36.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c37_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c37.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c38_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c38.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c39_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c39.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_c40_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_c40.c"
#include "mgga_xc_ml25_inc.c"

/* Integral 49-57: MN15 a terms (9 terms) */
#define XC_PREFIX           mn15_a01_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_a01.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_a02_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_a02.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_a03_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_a03.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_a04_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_a04.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_a05_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_a05.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_a06_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_a06.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_a07_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_a07.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_a08_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_a08.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_a09_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_a09.c"
#include "mgga_xc_ml25_inc.c"

/* Integral 58-66: MN15 b terms (9 terms) */
#define XC_PREFIX           mn15_b01_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_b01.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_b02_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_b02.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_b03_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_b03.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_b04_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_b04.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_b05_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_b05.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_b06_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_b06.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_b07_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_b07.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_b08_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_b08.c"
#include "mgga_xc_ml25_inc.c"

#define XC_PREFIX           mn15_b09_
#define XC_MAPLE2C_INC      "maple2c/mgga_exc/mgga_xc_ml25_mn15_b09.c"
#include "mgga_xc_ml25_inc.c"

#define N_INT 66

xc_mgga_funcs_variants *work_mgga_ml25[N_INT + 1] =
{
    &mn15_work_mgga,
    &lak_x_work_mgga,
    &lak_c_work_mgga,
    &lyp_c_t1_work_mgga,
    &lyp_c_t2_work_mgga,
    &lyp_c_t3_work_mgga,
    &lyp_c_t4_work_mgga,
    &lyp_c_t5_work_mgga,
    &lyp_c_t6_work_mgga,
    &mn15_c01_work_mgga,
    &mn15_c02_work_mgga,
    &mn15_c03_work_mgga,
    &mn15_c04_work_mgga,
    &mn15_c05_work_mgga,
    &mn15_c06_work_mgga,
    &mn15_c07_work_mgga,
    &mn15_c08_work_mgga,
    &mn15_c09_work_mgga,
    &mn15_c10_work_mgga,
    &mn15_c11_work_mgga,
    &mn15_c12_work_mgga,
    &mn15_c13_work_mgga,
    &mn15_c14_work_mgga,
    &mn15_c15_work_mgga,
    &mn15_c16_work_mgga,
    &mn15_c17_work_mgga,
    &mn15_c18_work_mgga,
    &mn15_c19_work_mgga,
    &mn15_c20_work_mgga,
    &mn15_c21_work_mgga,
    &mn15_c22_work_mgga,
    &mn15_c23_work_mgga,
    &mn15_c24_work_mgga,
    &mn15_c25_work_mgga,
    &mn15_c26_work_mgga,
    &mn15_c27_work_mgga,
    &mn15_c28_work_mgga,
    &mn15_c29_work_mgga,
    &mn15_c30_work_mgga,
    &mn15_c31_work_mgga,
    &mn15_c32_work_mgga,
    &mn15_c33_work_mgga,
    &mn15_c34_work_mgga,
    &mn15_c35_work_mgga,
    &mn15_c36_work_mgga,
    &mn15_c37_work_mgga,
    &mn15_c38_work_mgga,
    &mn15_c39_work_mgga,
    &mn15_c40_work_mgga,
    &mn15_a01_work_mgga,
    &mn15_a02_work_mgga,
    &mn15_a03_work_mgga,
    &mn15_a04_work_mgga,
    &mn15_a05_work_mgga,
    &mn15_a06_work_mgga,
    &mn15_a07_work_mgga,
    &mn15_a08_work_mgga,
    &mn15_a09_work_mgga,
    &mn15_b01_work_mgga,
    &mn15_b02_work_mgga,
    &mn15_b03_work_mgga,
    &mn15_b04_work_mgga,
    &mn15_b05_work_mgga,
    &mn15_b06_work_mgga,
    &mn15_b07_work_mgga,
    &mn15_b08_work_mgga,
    &mn15_b09_work_mgga,
};
