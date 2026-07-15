#include "ifxc_internal.h"

#include <string.h>

static ifxc_status
ifxc_get_dims_from_handle(const ifxc_func_type *func, ifxc_dimensions_t *dims)
{
  ifxc_status status = ifxc_func_dimensions(func, dims);
  return status;
}

static ifxc_status
ifxc_validate_input(const ifxc_dimensions_t *dims, const ifxc_input *input)
{
  if(dims == NULL || input == NULL){
    return IFXC_E_INVALID_ARGUMENT;
  }

  if(input->npoints > 0){
    if(input->rho == NULL || input->sigma == NULL || input->tau == NULL){
      return IFXC_E_MISSING_INPUT;
    }
  }

  return IFXC_OK;
}

static int
ifxc_var_component_count(const ifxc_dimensions_t *dims, ifxc_variable var)
{
  switch(var){
  case IFXC_VAR_RHO:
    return (int)dims->rho;
  case IFXC_VAR_SIGMA:
    return (int)dims->sigma;
  case IFXC_VAR_LAPL:
    return (int)dims->lapl;
  case IFXC_VAR_TAU:
    return (int)dims->tau;
  default:
    return -1;
  }
}

static ifxc_status
ifxc_validate_deriv_entry(
    const ifxc_dimensions_t *dims,
    unsigned int max_order,
    const ifxc_deriv_entry *entry,
    int require_output)
{
  unsigned int i;

  if(entry == NULL || (require_output && entry->out == NULL)){
    return IFXC_E_INVALID_ARGUMENT;
  }
  if(entry->order > max_order){
    return IFXC_E_UNSUPPORTED_DERIVATIVE;
  }

  if(entry->order == 0){
    return IFXC_OK;
  }
  if(entry->vars == NULL){
    return IFXC_E_INVALID_ARGUMENT;
  }

  for(i = 0; i < entry->order; ++i){
    if(entry->vars[i] > IFXC_VAR_TAU){
      return IFXC_E_INVALID_ARGUMENT;
    }
    if(i > 0 && entry->vars[i] < entry->vars[i - 1]){
      return IFXC_E_INVALID_ARGUMENT;
    }
    if(entry->vars[i] == IFXC_VAR_LAPL){
      return IFXC_E_UNSUPPORTED_VARIABLE;
    }
    if(ifxc_var_component_count(dims, entry->vars[i]) <= 0){
      return IFXC_E_UNSUPPORTED_VARIABLE;
    }
  }

  return IFXC_OK;
}

static ifxc_status
ifxc_entry_size(const ifxc_func_type *func, const ifxc_input *input,
                size_t nfeatures,
                unsigned int max_order,
                const ifxc_deriv_entry *entry,
                size_t *n_double)
{
  ifxc_dimensions_t dims;
  ifxc_status status;
  size_t count = 1;

  if(n_double == NULL){
    return IFXC_E_INVALID_ARGUMENT;
  }

  status = ifxc_get_dims_from_handle(func, &dims);
  if(status != IFXC_OK){
    return status;
  }
  status = ifxc_validate_input(&dims, input);
  if(status != IFXC_OK){
    return status;
  }
  status = ifxc_validate_deriv_entry(&dims, max_order, entry, 0);
  if(status != IFXC_OK){
    return status;
  }
  if(entry->target != IFXC_TARGET_LOCAL &&
     entry->target != IFXC_TARGET_INTEGRAL){
    return IFXC_E_INVALID_ARGUMENT;
  }

  if(entry->order == 0){
    if(entry->target == IFXC_TARGET_LOCAL){
      *n_double = nfeatures * input->npoints;
      return IFXC_OK;
    }
    if(entry->target == IFXC_TARGET_INTEGRAL){
      *n_double = nfeatures;
      return IFXC_OK;
    }
    return IFXC_E_INVALID_ARGUMENT;
  }

  status = ifxc_derivative_component_count(&dims, entry, &count);
  if(status != IFXC_OK){
    return status;
  }

  *n_double = nfeatures * count * input->npoints;
  return IFXC_OK;
}

ifxc_status
ifxc_output_size(const ifxc_func_type *func, const ifxc_input *input,
                 const ifxc_deriv_entry *entry, size_t *n_double)
{
  const ifxc_handle_impl *impl = (const ifxc_handle_impl *)((func != NULL) ? func->impl : NULL);
  if(impl == NULL){
    return IFXC_E_NOT_INITIALIZED;
  }

  return ifxc_entry_size(func, input, impl->nfeatures,
                         impl->max_deriv_order, entry, n_double);
}

static ifxc_status
ifxc_validate_eval_request(const ifxc_func_type *func, const ifxc_input *input,
                           unsigned int max_order,
                           size_t nentries, const ifxc_deriv_entry *entries,
                           ifxc_dimensions_t *dims)
{
  size_t i;
  ifxc_status status;

  if(func == NULL || input == NULL || dims == NULL){
    return IFXC_E_INVALID_ARGUMENT;
  }

  status = ifxc_get_dims_from_handle(func, dims);
  if(status != IFXC_OK){
    return status;
  }

  status = ifxc_validate_input(dims, input);
  if(status != IFXC_OK){
    return status;
  }

  if(nentries == 0 || entries == NULL){
    return IFXC_E_INVALID_ARGUMENT;
  }

  for(i = 0; i < nentries; ++i){
    status = ifxc_validate_deriv_entry(dims, max_order, &entries[i], 1);
    if(status != IFXC_OK){
      return status;
    }
    if(entries[i].target == IFXC_TARGET_INTEGRAL && input->weights == NULL){
      return IFXC_E_MISSING_INPUT;
    }
    if(entries[i].target != IFXC_TARGET_LOCAL &&
       entries[i].target != IFXC_TARGET_INTEGRAL){
      return IFXC_E_INVALID_ARGUMENT;
    }
  }

  return IFXC_OK;
}

ifxc_status
ifxc_eval(const ifxc_func_type *func, const ifxc_input *input,
          size_t nentries, const ifxc_deriv_entry *entries)
{
  ifxc_dimensions_t dims;
  ifxc_status status;
  const ifxc_handle_impl *impl;

  impl = (const ifxc_handle_impl *)((func != NULL) ? func->impl : NULL);
  if(impl == NULL){
    return IFXC_E_NOT_INITIALIZED;
  }

  status = ifxc_validate_eval_request(func, input, impl->max_deriv_order, nentries, entries, &dims);
  if(status != IFXC_OK){
    return status;
  }

  (void)dims;
  (void)impl;
  return ifxc_ml25_eval(func, input, nentries, entries);
}

ifxc_status
ifxc_eval_first_derivatives_contracted(
    const ifxc_func_type *func,
    const ifxc_input *input,
    const double *coeffs,
    double *d_rho,
    double *d_sigma,
    double *d_tau)
{
  ifxc_dimensions_t dims;
  ifxc_status status;
  const ifxc_handle_impl *impl;

  impl = (const ifxc_handle_impl *)((func != NULL) ? func->impl : NULL);
  if(impl == NULL){
    return IFXC_E_NOT_INITIALIZED;
  }
  if(input == NULL || coeffs == NULL || d_rho == NULL ||
     d_sigma == NULL || d_tau == NULL){
    return IFXC_E_INVALID_ARGUMENT;
  }
  status = ifxc_get_dims_from_handle(func, &dims);
  if(status != IFXC_OK){
    return status;
  }
  status = ifxc_validate_input(&dims, input);
  if(status != IFXC_OK){
    return status;
  }
  switch(impl->feature_set){
  case IFXC_FEATURE_SET_ML25:
    return ifxc_ml25_eval_first_derivatives_contracted(
        impl, input, coeffs, d_rho, d_sigma, d_tau);
  default:
    return IFXC_E_UNKNOWN_FEATURE_SET;
  }
}

ifxc_status
ifxc_eval_ml25_first_derivatives_contracted(
    const ifxc_func_type *func,
    const ifxc_input *input,
    const double *coeffs,
    double *d_rho,
    double *d_sigma,
    double *d_tau)
{
  return ifxc_eval_first_derivatives_contracted(
      func, input, coeffs, d_rho, d_sigma, d_tau);
}
