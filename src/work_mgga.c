/*
 Copyright (C) 2006-2018 M.A.L. Marques
 Copyright (C) 2019 X. Andrade

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/**
 * @file work_mgga.c
 * @brief This file is to be included in MGGA functionals.
 */

#ifdef XC_PREFIX
  #ifndef ADD_PREFIX
  #define ADD_PREFIX_(a, b) a ## b
  #define ADD_PREFIX(a, b) ADD_PREFIX_(a, b)
  #endif
#endif

/* define auxiliary functions to NULL in case they are not available */
#if defined(XC_DONT_COMPILE_EXC) || maple2c_order < 0 || defined(XC_NO_EXC)
#define work_mgga_exc_unpol NULL
#define work_mgga_exc_pol NULL
#else
  #ifdef XC_PREFIX
  #define work_mgga_exc_unpol ADD_PREFIX(XC_PREFIX, work_mgga_exc_unpol)
  #define work_mgga_exc_pol   ADD_PREFIX(XC_PREFIX, work_mgga_exc_pol)
  #endif
#define ORDER_TXT exc
#define SPIN_TXT  unpol
#include "work_mgga_inc.c"
#undef SPIN_TXT
#define SPIN_TXT  pol
#include "work_mgga_inc.c"
#undef SPIN_TXT
#undef ORDER_TXT
#endif

#if defined(XC_DONT_COMPILE_VXC) || maple2c_order < 1
#define work_mgga_vxc_unpol NULL
#define work_mgga_vxc_pol NULL
#else
  #ifdef XC_PREFIX
  #define work_mgga_vxc_unpol ADD_PREFIX(XC_PREFIX, work_mgga_vxc_unpol)
  #define work_mgga_vxc_pol   ADD_PREFIX(XC_PREFIX, work_mgga_vxc_pol)
  #endif
#define ORDER_TXT vxc
#define SPIN_TXT  unpol
#include "work_mgga_inc.c"
#undef SPIN_TXT
#define SPIN_TXT  pol
#include "work_mgga_inc.c"
#undef SPIN_TXT
#undef ORDER_TXT
#endif

#if defined(XC_DONT_COMPILE_FXC) || maple2c_order < 2
#define work_mgga_fxc_unpol NULL
#define work_mgga_fxc_pol NULL
#else
  #ifdef XC_PREFIX
  #define work_mgga_fxc_unpol ADD_PREFIX(XC_PREFIX, work_mgga_fxc_unpol)
  #define work_mgga_fxc_pol   ADD_PREFIX(XC_PREFIX, work_mgga_fxc_pol)
  #endif
#define ORDER_TXT fxc
#define SPIN_TXT  unpol
#include "work_mgga_inc.c"
#undef SPIN_TXT
#define SPIN_TXT  pol
#include "work_mgga_inc.c"
#undef SPIN_TXT
#undef ORDER_TXT
#endif

#if defined(XC_DONT_COMPILE_KXC) || maple2c_order < 3
#define work_mgga_kxc_unpol NULL
#define work_mgga_kxc_pol NULL
#else
  #ifdef XC_PREFIX
  #define work_mgga_kxc_unpol ADD_PREFIX(XC_PREFIX, work_mgga_kxc_unpol)
  #define work_mgga_kxc_pol   ADD_PREFIX(XC_PREFIX, work_mgga_kxc_pol)
  #endif
#define ORDER_TXT kxc
#define SPIN_TXT  unpol
#include "work_mgga_inc.c"
#undef SPIN_TXT
#define SPIN_TXT  pol
#include "work_mgga_inc.c"
#undef SPIN_TXT
#undef ORDER_TXT
#endif

#if defined(XC_DONT_COMPILE_LXC) || maple2c_order < 4
#define work_mgga_lxc_unpol NULL
#define work_mgga_lxc_pol NULL
#else
  #ifdef XC_PREFIX
  #define work_mgga_lxc_unpol ADD_PREFIX(XC_PREFIX, work_mgga_lxc_unpol)
  #define work_mgga_lxc_pol   ADD_PREFIX(XC_PREFIX, work_mgga_lxc_pol)
  #endif
#define ORDER_TXT lxc
#define SPIN_TXT  unpol
#include "work_mgga_inc.c"
#undef SPIN_TXT
#define SPIN_TXT  pol
#include "work_mgga_inc.c"
#undef SPIN_TXT
#undef ORDER_TXT
#endif

#ifdef XC_PREFIX
#define work_mgga ADD_PREFIX(XC_PREFIX, work_mgga)
#endif

/* we construct a structure containing all variants */
static xc_mgga_funcs_variants work_mgga =
  {
   {work_mgga_exc_unpol, work_mgga_vxc_unpol, work_mgga_fxc_unpol, work_mgga_kxc_unpol, work_mgga_lxc_unpol},
   {work_mgga_exc_pol,   work_mgga_vxc_pol,   work_mgga_fxc_pol,   work_mgga_kxc_pol,   work_mgga_lxc_pol}
  };

#ifdef XC_PREFIX
#undef work_mgga
#undef work_mgga_exc_unpol
#undef work_mgga_vxc_unpol
#undef work_mgga_fxc_unpol
#undef work_mgga_kxc_unpol
#undef work_mgga_lxc_unpol
#undef work_mgga_exc_pol
#undef work_mgga_vxc_pol
#undef work_mgga_fxc_pol
#undef work_mgga_kxc_pol
#undef work_mgga_lxc_pol
#endif
