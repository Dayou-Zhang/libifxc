/* SPDX-License-Identifier: MPL-2.0 */

#ifndef IFXC_MATH_H
#define IFXC_MATH_H

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ifxc_mgga.h"

#if defined(__GNUC__) || defined(__clang__)
#define FORCE_INLINE __attribute__((always_inline)) inline
#else
#define FORCE_INLINE inline
#endif

#define GPU_DEVICE_FUNCTION

#define IFXC_MGGA_FLAGS_I_HAVE_EXC IFXC_MGGA_FLAGS_HAVE_EXC
#define IFXC_MGGA_FLAGS_I_HAVE_VXC IFXC_MGGA_FLAGS_HAVE_VXC
#define IFXC_MGGA_FLAGS_I_HAVE_FXC IFXC_MGGA_FLAGS_HAVE_FXC
#define IFXC_MGGA_FLAGS_I_HAVE_KXC IFXC_MGGA_FLAGS_HAVE_KXC
#define IFXC_MGGA_FLAGS_I_HAVE_LXC IFXC_MGGA_FLAGS_HAVE_LXC

#define IFXC_MGGA_FLAGS_I_HAVE_ALL \
  (IFXC_MGGA_FLAGS_I_HAVE_EXC | IFXC_MGGA_FLAGS_I_HAVE_VXC | \
   IFXC_MGGA_FLAGS_I_HAVE_FXC | IFXC_MGGA_FLAGS_I_HAVE_KXC | \
   IFXC_MGGA_FLAGS_I_HAVE_LXC)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif

#define POW_1_2(x) sqrt(x)
#define POW_1_4(x) sqrt(sqrt(x))
#define POW_2(x) ((x) * (x))
#define POW_3_2(x) ((x) * sqrt(x))

#if defined(__GNUC__) || defined(__clang__)
#define POW_1_3(x) cbrt(x)
#else
#define POW_1_3(x) pow((x), 1.0 / 3.0)
#endif

#define my_piecewise3(c, x1, x2) ((c) ? (x1) : (x2))
#define my_piecewise5(c1, x1, c2, x2, x3) ((c1) ? (x1) : ((c2) ? (x2) : (x3)))

#define M_SQRTPI 1.772453850905516027298167483341145182798L
#define M_CBRTPI 1.464591887561523263020142527263790391739L
#define M_SQRT3 1.732050807568877293527446341505872366943L
#define M_CBRT2 1.259921049894873164767210607278350570L
#define M_CBRT3 1.442249570307408382321638310780109588392L
#define M_CBRT4 1.587401051968199474751705639272308260391L
#define M_CBRT5 1.709975946676696989353108872543860109868L
#define M_CBRT6 1.817120592832139658891211756327260502428L
#define M_CBRT7 1.912931182772389101199116839548760282862L
#define M_CBRT9 2.080083823051904114530056824357885386338L

#ifndef m_min
#define m_min(x, y) (((x) < (y)) ? (x) : (y))
#endif
#ifndef m_max
#define m_max(x, y) (((x) < (y)) ? (y) : (x))
#endif

#endif /* IFXC_MATH_H */
