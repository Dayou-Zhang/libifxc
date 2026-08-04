/* SPDX-License-Identifier: MPL-2.0 */

#include "ifxc_internal.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "ifxc_mgga.h"
#include "ifxc_vector_mgga.h"

ifxc_mgga_funcs_variants *ifxc_ml26_combined_work_mgga(void);
int ifxc_ml26_has_combined_order(unsigned int order);
ifxc_status ifxc_ml26_combined_eval_contracted_component_major(
    const ifxc_mgga_func_type *p,
    size_t npoints,
    const double *rho,
    const double *sigma,
    const double *tau,
    double *d_rho,
    double *d_sigma,
    double *d_tau);
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
           IFXC_MGGA_FLAGS_HAVE_ALL,
  .dens_threshold = 1e-15
};

static void
ifxc_ml26_prepare_private_func(ifxc_mgga_func_type *p, ifxc_nspin nspin)
{
  memset(p, 0, sizeof(*p));
  p->info = &IFXC_ML26_PRIVATE_INFO;
  p->nspin = (nspin == IFXC_POLARIZED)
      ? IFXC_MGGA_POLARIZED
      : IFXC_MGGA_UNPOLARIZED;
  p->dens_threshold = p->info->dens_threshold;
  p->sigma_threshold = pow(p->dens_threshold, 4.0 / 3.0);
  p->zeta_threshold = DBL_EPSILON;
  p->tau_threshold = 1e-20;
}

static ifxc_status
ifxc_ml26_replace_ml25_derivative_prefix(
    const ifxc_handle_impl *impl,
    const ifxc_input *input,
    size_t nentries,
    const ifxc_deriv_entry *entries)
{
  ifxc_handle_impl prefix_impl = *impl;
  ifxc_func_type prefix_func;
  size_t entry_index;

  prefix_impl.feature_set = IFXC_FEATURE_SET_ML25;
  prefix_impl.nfeatures = IFXC_ML25_NFEATURES;
  prefix_impl.max_deriv_order = IFXC_ML25_GENERATED_MAX_ORDER;
  prefix_func.impl = &prefix_impl;

  for(entry_index = 0; entry_index < nentries; ++entry_index){
    const ifxc_deriv_entry *entry = &entries[entry_index];
    ifxc_deriv_entry prefix_entry;
    size_t ncomponents;
    size_t nvalues;
    double *prefix;
    size_t point;
    size_t feature;
    ifxc_status status;

    if(entry->order == 0){
      continue;
    }
    status = ifxc_derivative_component_count(
        &impl->dims, entry, &ncomponents);
    if(status != IFXC_OK){
      return status;
    }
    nvalues = input->npoints * IFXC_ML25_NFEATURES * ncomponents;
    prefix = (double *)calloc(nvalues, sizeof(double));
    if(prefix == NULL && nvalues > 0){
      return IFXC_E_ALLOCATION;
    }
    prefix_entry = *entry;
    prefix_entry.out = prefix;
    status = ifxc_ml25_eval(
        &prefix_func, input, 1, &prefix_entry);
    if(status != IFXC_OK){
      free(prefix);
      return status;
    }
    for(point = 0; point < input->npoints; ++point){
      for(feature = 0; feature < IFXC_ML25_NFEATURES; ++feature){
        size_t component;
        for(component = 0; component < ncomponents; ++component){
          size_t src =
              (point * IFXC_ML25_NFEATURES + feature) * ncomponents
              + component;
          size_t dst =
              (point * IFXC_ML26_NFEATURES + feature) * ncomponents
              + component;
          entry->out[dst] = prefix[src];
        }
      }
    }
    free(prefix);
  }
  return IFXC_OK;
}

static ifxc_status
ifxc_ml26_eval_combined(
    const ifxc_handle_impl *impl,
    const ifxc_input *input,
    size_t nentries,
    const ifxc_deriv_entry *entries,
    unsigned int max_order)
{
  ifxc_mgga_dimensions private_dims;
  double *rho_tm = NULL;
  double *sigma_tm = NULL;
  double *lapl_tm = NULL;
  double *tau_tm = NULL;
  ifxc_mgga_out_params out = {0};
  ifxc_mgga_func_type p = {0};
  ifxc_status status;
  size_t entry_index;
  ifxc_mgga_funcs selected_work;

  ifxc_vector_mgga_set_private_dimensions(
      (impl->nspin == IFXC_UNPOLARIZED)
          ? IFXC_MGGA_UNPOLARIZED
          : IFXC_MGGA_POLARIZED,
      &private_dims);
  ifxc_vector_mgga_zero_dimensions_above_order(&private_dims, max_order);
  ifxc_vector_mgga_scale_output_dimensions(&private_dims, impl->nfeatures);

  status = ifxc_vector_mgga_transpose_inputs(
      &impl->dims, input, &rho_tm, &sigma_tm, &lapl_tm, &tau_tm);
  if(status != IFXC_OK){
    goto cleanup;
  }
  status = ifxc_vector_mgga_allocate_private_outputs(
      &private_dims, input->npoints, &out);
  if(status != IFXC_OK){
    goto cleanup;
  }

  ifxc_ml26_prepare_private_func(&p, impl->nspin);
  p.dim = private_dims;
  selected_work = (impl->nspin == IFXC_UNPOLARIZED)
      ? ifxc_ml26_combined_work_mgga()->unpol[max_order]
      : ifxc_ml26_combined_work_mgga()->pol[max_order];
  if(selected_work == NULL){
    status = IFXC_E_INTERNAL;
    goto cleanup;
  }
  selected_work(
      &p, input->npoints, rho_tm, sigma_tm, lapl_tm, tau_tm, &out);

  for(entry_index = 0; entry_index < nentries; ++entry_index){
    const ifxc_deriv_entry *entry = &entries[entry_index];
    const double *src =
        ifxc_vector_mgga_select_source_buffer(&out, entry);
    if(src == NULL){
      status = IFXC_E_INTERNAL;
      goto cleanup;
    }
    status = ifxc_vector_mgga_copy_combined_entry(
        input, &impl->dims, impl->nfeatures, entry, src, 1.0);
    if(status != IFXC_OK){
      goto cleanup;
    }
  }
  status = ifxc_ml26_replace_ml25_derivative_prefix(
      impl, input, nentries, entries);

cleanup:
  ifxc_vector_mgga_free_private_outputs(&out);
  ifxc_vector_mgga_free_inputs(rho_tm, sigma_tm, lapl_tm, tau_tm);
  return status;
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
  unsigned int max_order;

  if(impl == NULL || input == NULL || entries == NULL){
    return IFXC_E_INVALID_ARGUMENT;
  }
  if(impl->feature_set != IFXC_FEATURE_SET_ML26){
    return IFXC_E_UNKNOWN_FEATURE_SET;
  }
  if(nentries == 1 && entries[0].target == IFXC_TARGET_INTEGRAL &&
     entries[0].order == 0){
    return ifxc_ml26_eval_integral_order0(impl, input, entries[0].out);
  }

  max_order = ifxc_vector_mgga_max_requested_order(nentries, entries);
  if(max_order <= IFXC_ML26_GENERATED_MAX_ORDER &&
     ifxc_ml26_has_combined_order(max_order)){
    return ifxc_ml26_eval_combined(
        impl, input, nentries, entries, max_order);
  }
  return IFXC_E_UNSUPPORTED_DERIVATIVE;
}

ifxc_status
ifxc_ml26_eval_first_derivatives_contracted(
    const ifxc_handle_impl *impl,
    const ifxc_input *input,
    const double *coeffs,
    double *d_rho,
    double *d_sigma,
    double *d_tau)
{
  ifxc_mgga_dimensions private_dims;
  ifxc_mgga_func_type p = {0};

  if(impl == NULL || input == NULL || coeffs == NULL ||
     d_rho == NULL || d_sigma == NULL || d_tau == NULL){
    return IFXC_E_INVALID_ARGUMENT;
  }
  if(impl->feature_set != IFXC_FEATURE_SET_ML26){
    return IFXC_E_UNKNOWN_FEATURE_SET;
  }
  if(!ifxc_ml26_has_combined_order(1)){
    return IFXC_E_UNSUPPORTED_DERIVATIVE;
  }

  ifxc_vector_mgga_set_private_dimensions(
      (impl->nspin == IFXC_UNPOLARIZED)
          ? IFXC_MGGA_UNPOLARIZED
          : IFXC_MGGA_POLARIZED,
      &private_dims);
  ifxc_vector_mgga_zero_dimensions_above_order(&private_dims, 1);
  private_dims.zk = 0;
  ifxc_ml26_prepare_private_func(&p, impl->nspin);
  p.dim = private_dims;
  p.feature_coeffs = coeffs;

  return ifxc_ml26_combined_eval_contracted_component_major(
      &p, input->npoints, input->rho, input->sigma, input->tau,
      d_rho, d_sigma, d_tau);
}

ifxc_status
ifxc_ml26_eval_integral_order0(
    const ifxc_handle_impl *impl,
    const ifxc_input *input,
    double *features)
{
  ifxc_mgga_dimensions private_dims;
  ifxc_mgga_func_type p = {0};

  if(impl == NULL || input == NULL || features == NULL ||
     input->weights == NULL){
    return IFXC_E_INVALID_ARGUMENT;
  }
  if(impl->feature_set != IFXC_FEATURE_SET_ML26){
    return IFXC_E_UNKNOWN_FEATURE_SET;
  }

  ifxc_vector_mgga_set_private_dimensions(
      (impl->nspin == IFXC_UNPOLARIZED)
          ? IFXC_MGGA_UNPOLARIZED
          : IFXC_MGGA_POLARIZED,
      &private_dims);
  ifxc_vector_mgga_zero_dimensions_above_order(&private_dims, 0);
  ifxc_vector_mgga_scale_output_dimensions(
      &private_dims, impl->nfeatures);
  ifxc_ml26_prepare_private_func(&p, impl->nspin);
  p.dim = private_dims;

  return ifxc_ml26_combined_eval_order0_component_major(
      &p, input->npoints, input->rho, input->sigma, input->tau,
      input->weights, NULL, features);
}
