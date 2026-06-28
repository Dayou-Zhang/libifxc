#include "util.h"

#define XC_PREFIX ml25_combined_
#define XC_MAPLE2C_INC "maple2c/if_mgga/mgga_xc_ml25.c"
#include "mgga_xc_ml25_inc.c"

ifxc_mgga_funcs_variants *
ifxc_ml25_combined_work_mgga(void)
{
  return &ml25_combined_work_mgga;
}

int
ifxc_ml25_has_combined_order(unsigned int order)
{
  if(order >= 5){
    return 0;
  }
  return ml25_combined_work_mgga.unpol[order] != NULL &&
         ml25_combined_work_mgga.pol[order] != NULL;
}
