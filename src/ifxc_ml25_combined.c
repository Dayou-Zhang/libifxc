#include "util.h"

#define XC_PREFIX ml25_combined_
#define XC_MAPLE2C_INC "maple2c/mgga_exc/mgga_xc_ml25.c"
#include "mgga_xc_ml25_inc.c"

xc_mgga_funcs_variants *
ifxc_ml25_combined_work_mgga(void)
{
  return &ml25_combined_work_mgga;
}

int
ifxc_ml25_has_combined_order0(void)
{
  return ml25_combined_work_mgga.unpol[0] != NULL &&
         ml25_combined_work_mgga.pol[0] != NULL;
}
