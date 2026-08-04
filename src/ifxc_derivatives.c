/* SPDX-License-Identifier: MPL-2.0 */

#include "ifxc_internal.h"

static size_t
ifxc_var_dimension(const ifxc_dimensions_t *dims, ifxc_variable var)
{
  switch(var){
  case IFXC_VAR_RHO:
    return dims->rho;
  case IFXC_VAR_SIGMA:
    return dims->sigma;
  case IFXC_VAR_LAPL:
    return dims->lapl;
  case IFXC_VAR_TAU:
    return dims->tau;
  default:
    return 0;
  }
}

static size_t
ifxc_binomial(size_t n, size_t k)
{
  size_t i;
  size_t result = 1;

  if(k > n){
    return 0;
  }
  if(k > n - k){
    k = n - k;
  }

  for(i = 1; i <= k; ++i){
    result = (result * (n - k + i)) / i;
  }

  return result;
}

ifxc_status
ifxc_derivative_component_count(
    const ifxc_dimensions_t *dims,
    const ifxc_deriv_entry *entry,
    size_t *ncomponents)
{
  unsigned int counts[4] = {0, 0, 0, 0};
  size_t total = 1;
  unsigned int i;

  if(dims == NULL || entry == NULL || ncomponents == NULL){
    return IFXC_E_INVALID_ARGUMENT;
  }

  if(entry->order == 0){
    *ncomponents = 1;
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
    counts[entry->vars[i]] += 1;
  }

  for(i = 0; i <= (unsigned int)IFXC_VAR_TAU; ++i){
    size_t dim;

    if(counts[i] == 0){
      continue;
    }

    dim = ifxc_var_dimension(dims, (ifxc_variable)i);
    if(dim == 0){
      return IFXC_E_UNSUPPORTED_VARIABLE;
    }
    total *= ifxc_binomial(dim + (size_t)counts[i] - 1, counts[i]);
  }

  *ncomponents = total;
  return IFXC_OK;
}
