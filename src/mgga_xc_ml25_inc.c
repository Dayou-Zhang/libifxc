#ifndef ADD_PREFIX
#define ADD_PREFIX_(a, b) a ## b
#define ADD_PREFIX(a, b) ADD_PREFIX_(a, b)
#endif

#define func_exc_unpol      ADD_PREFIX(XC_PREFIX, func_exc_unpol)
#define func_vxc_unpol      ADD_PREFIX(XC_PREFIX, func_vxc_unpol)
#define func_fxc_unpol      ADD_PREFIX(XC_PREFIX, func_fxc_unpol)
#define func_kxc_unpol      ADD_PREFIX(XC_PREFIX, func_kxc_unpol)
#define func_lxc_unpol      ADD_PREFIX(XC_PREFIX, func_lxc_unpol)
#define func_exc_pol        ADD_PREFIX(XC_PREFIX, func_exc_pol)
#define func_vxc_pol        ADD_PREFIX(XC_PREFIX, func_vxc_pol)
#define func_fxc_pol        ADD_PREFIX(XC_PREFIX, func_fxc_pol)
#define func_kxc_pol        ADD_PREFIX(XC_PREFIX, func_kxc_pol)
#define func_lxc_pol        ADD_PREFIX(XC_PREFIX, func_lxc_pol)
#include XC_MAPLE2C_INC
#include "ifxc_work_mgga.c"
#undef func_exc_unpol
#undef func_vxc_unpol
#undef func_fxc_unpol
#undef func_kxc_unpol
#undef func_lxc_unpol
#undef func_exc_pol
#undef func_vxc_pol
#undef func_fxc_pol
#undef func_kxc_pol
#undef func_lxc_pol
#undef XC_PREFIX
#undef XC_MAPLE2C_INC
