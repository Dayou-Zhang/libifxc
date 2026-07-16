#include "ifxc_internal.h"

#include <float.h>
#include <math.h>
#include <string.h>

#include "ifxc_mgga.h"

ifxc_status ifxc_ml26_combined_eval_order0_component_major(
    const ifxc_mgga_func_type *p,
    size_t npoints,
    const double *rho,
    const double *sigma,
    const double *tau,
    const double *weights,
    double *local_values,
    double *features);

static const ifxc_mgga_func_info IFXC_ML26_PRIVATE_INFO = {
  .flags = IFXC_MGGA_FLAGS_3D | IFXC_MGGA_FLAGS_NEEDS_TAU |
           IFXC_MGGA_FLAGS_HAVE_EXC,
  .dens_threshold = 1e-15
};

static void
ifxc_ml26_prepare_private_func(
    ifxc_mgga_func_type *p,
    ifxc_nspin nspin)
{
  memset(p, 0, sizeof(*p));
  p->info = &IFXC_ML26_PRIVATE_INFO;
  p->nspin = (nspin == IFXC_POLARIZED)
      ? IFXC_MGGA_POLARIZED
      : IFXC_MGGA_UNPOLARIZED;
  p->dim.zk = IFXC_ML26_NFEATURES;
  if(p->nspin == IFXC_MGGA_POLARIZED){
    p->dim.rho = 2;
    p->dim.sigma = 3;
    p->dim.tau = 2;
  }else{
    p->dim.rho = 1;
    p->dim.sigma = 1;
    p->dim.tau = 1;
  }
  p->dens_threshold = p->info->dens_threshold;
  p->sigma_threshold = pow(p->dens_threshold, 4.0 / 3.0);
  p->zeta_threshold = DBL_EPSILON;
  p->tau_threshold = 1e-20;
}

ifxc_status
ifxc_ml26_eval(
    const ifxc_func_type *func,
    const ifxc_input *input,
    size_t nentries,
    const ifxc_deriv_entry *entries)
{
  const ifxc_handle_impl *impl =
      (const ifxc_handle_impl *)((func != NULL) ? func->impl : NULL);
  ifxc_mgga_func_type p;
  size_t entry_index;

  if(impl == NULL || input == NULL || entries == NULL){
    return IFXC_E_INVALID_ARGUMENT;
  }
  if(impl->feature_set != IFXC_FEATURE_SET_ML26){
    return IFXC_E_UNKNOWN_FEATURE_SET;
  }

  ifxc_ml26_prepare_private_func(&p, impl->nspin);
  for(entry_index = 0; entry_index < nentries; ++entry_index){
    const ifxc_deriv_entry *entry = &entries[entry_index];
    if(entry->order != 0){
      return IFXC_E_UNSUPPORTED_DERIVATIVE;
    }
    if(entry->target == IFXC_TARGET_LOCAL){
      ifxc_status status = ifxc_ml26_combined_eval_order0_component_major(
          &p, input->npoints, input->rho, input->sigma, input->tau,
          NULL, entry->out, NULL);
      if(status != IFXC_OK){
        return status;
      }
    }else if(entry->target == IFXC_TARGET_INTEGRAL){
      ifxc_status status = ifxc_ml26_combined_eval_order0_component_major(
          &p, input->npoints, input->rho, input->sigma, input->tau,
          input->weights, NULL, entry->out);
      if(status != IFXC_OK){
        return status;
      }
    }else{
      return IFXC_E_INVALID_ARGUMENT;
    }
  }

  return IFXC_OK;
}
