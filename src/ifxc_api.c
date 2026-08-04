/* SPDX-License-Identifier: MPL-2.0 */

#include "ifxc_internal.h"

#include <stdlib.h>
#include <string.h>

static const char IFXC_VERSION_STRING[] = "0.4.0";

const char *
ifxc_version_string(void)
{
  return IFXC_VERSION_STRING;
}

int
ifxc_version_major(void)
{
  return IFXC_VERSION_MAJOR;
}

int
ifxc_version_minor(void)
{
  return IFXC_VERSION_MINOR;
}

int
ifxc_version_patch(void)
{
  return IFXC_VERSION_PATCH;
}

static ifxc_status
ifxc_validate_feature_set(int feature_set)
{
  switch(feature_set){
  case IFXC_FEATURE_SET_ML25:
  case IFXC_FEATURE_SET_ML26:
    return IFXC_OK;
  default:
    return IFXC_E_UNKNOWN_FEATURE_SET;
  }
}

static ifxc_status
ifxc_validate_nspin(ifxc_nspin nspin)
{
  switch(nspin){
  case IFXC_UNPOLARIZED:
  case IFXC_POLARIZED:
    return IFXC_OK;
  default:
    return IFXC_E_INVALID_NSPIN;
  }
}

ifxc_status
ifxc_init(ifxc_func_type *func, int feature_set, ifxc_nspin nspin)
{
  if(func == NULL){
    return IFXC_E_INVALID_ARGUMENT;
  }

  memset(func, 0, sizeof(*func));

  if(ifxc_validate_feature_set(feature_set) != IFXC_OK){
    return IFXC_E_UNKNOWN_FEATURE_SET;
  }
  if(ifxc_validate_nspin(nspin) != IFXC_OK){
    return IFXC_E_INVALID_NSPIN;
  }

  ifxc_handle_impl *impl = (ifxc_handle_impl *)calloc(1, sizeof(*impl));
  if(impl == NULL){
    return IFXC_E_ALLOCATION;
  }

  impl->feature_set = feature_set;
  impl->nspin = nspin;
  switch(feature_set){
  case IFXC_FEATURE_SET_ML25:
    impl->nfeatures = IFXC_ML25_NFEATURES;
    impl->max_deriv_order = IFXC_ML25_GENERATED_MAX_ORDER;
    break;
  case IFXC_FEATURE_SET_ML26:
    impl->nfeatures = IFXC_ML26_NFEATURES;
    impl->max_deriv_order = IFXC_ML26_GENERATED_MAX_ORDER;
    break;
  default:
    free(impl);
    return IFXC_E_UNKNOWN_FEATURE_SET;
  }

  ifxc_status status = ifxc_dimensions(feature_set, nspin, &impl->dims);
  if(status != IFXC_OK){
    free(impl);
    return status;
  }

  func->impl = impl;
  return IFXC_OK;
}

void
ifxc_end(ifxc_func_type *func)
{
  if(func == NULL){
    return;
  }

  free(func->impl);
  func->impl = NULL;
}

static const ifxc_handle_impl *
ifxc_get_impl(const ifxc_func_type *func)
{
  if(func == NULL){
    return NULL;
  }
  return (const ifxc_handle_impl *)func->impl;
}

ifxc_status
ifxc_nfeatures(const ifxc_func_type *func, size_t *nfeatures)
{
  const ifxc_handle_impl *impl = ifxc_get_impl(func);
  if(impl == NULL){
    return IFXC_E_NOT_INITIALIZED;
  }
  if(nfeatures == NULL){
    return IFXC_E_INVALID_ARGUMENT;
  }

  *nfeatures = impl->nfeatures;
  return IFXC_OK;
}

ifxc_status
ifxc_max_deriv_order(const ifxc_func_type *func, unsigned int *max_order)
{
  const ifxc_handle_impl *impl = ifxc_get_impl(func);
  if(impl == NULL){
    return IFXC_E_NOT_INITIALIZED;
  }
  if(max_order == NULL){
    return IFXC_E_INVALID_ARGUMENT;
  }

  *max_order = impl->max_deriv_order;
  return IFXC_OK;
}

ifxc_status
ifxc_func_dimensions(const ifxc_func_type *func, ifxc_dimensions_t *dims)
{
  const ifxc_handle_impl *impl = ifxc_get_impl(func);
  if(impl == NULL){
    return IFXC_E_NOT_INITIALIZED;
  }
  if(dims == NULL){
    return IFXC_E_INVALID_ARGUMENT;
  }

  *dims = impl->dims;
  return IFXC_OK;
}
