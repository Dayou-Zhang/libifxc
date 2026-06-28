#ifndef IFXC_PRIVATE_XC_H
#define IFXC_PRIVATE_XC_H

/*
 * Temporary compatibility names for retained generated code. New libifxc
 * backend code should include ifxc_mgga.h and use ifxc_ private types.
 */

#include "ifxc_mgga.h"

#define XC_UNPOLARIZED IFXC_MGGA_UNPOLARIZED
#define XC_POLARIZED IFXC_MGGA_POLARIZED

#define XC_FLAGS_HAVE_EXC IFXC_MGGA_FLAGS_HAVE_EXC
#define XC_FLAGS_HAVE_VXC IFXC_MGGA_FLAGS_HAVE_VXC
#define XC_FLAGS_HAVE_FXC IFXC_MGGA_FLAGS_HAVE_FXC
#define XC_FLAGS_HAVE_KXC IFXC_MGGA_FLAGS_HAVE_KXC
#define XC_FLAGS_HAVE_LXC IFXC_MGGA_FLAGS_HAVE_LXC
#define XC_FLAGS_3D IFXC_MGGA_FLAGS_3D
#define XC_FLAGS_NEEDS_LAPLACIAN IFXC_MGGA_FLAGS_NEEDS_LAPLACIAN
#define XC_FLAGS_NEEDS_TAU IFXC_MGGA_FLAGS_NEEDS_TAU
#define XC_FLAGS_ENFORCE_FHC IFXC_MGGA_FLAGS_ENFORCE_FHC

#define XC_FLAGS_HAVE_ALL IFXC_MGGA_FLAGS_HAVE_ALL

typedef ifxc_mgga_func_info xc_func_info_type;
typedef ifxc_mgga_dimensions xc_dimensions;
typedef ifxc_mgga_out_params xc_mgga_out_params;
typedef ifxc_mgga_funcs xc_mgga_funcs;
typedef ifxc_mgga_funcs_variants xc_mgga_funcs_variants;
typedef ifxc_mgga_func_type xc_func_type;

#define XC(func) xc_##func

#endif /* IFXC_PRIVATE_XC_H */
