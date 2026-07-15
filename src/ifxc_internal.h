#ifndef IFXC_INTERNAL_H
#define IFXC_INTERNAL_H

#include "ifxc.h"

#define IFXC_ML25_GENERATED_MAX_ORDER 3u

typedef struct ifxc_handle_impl {
  int feature_set;
  ifxc_nspin nspin;
  ifxc_dimensions_t dims;
  size_t nfeatures;
  unsigned int max_deriv_order;
} ifxc_handle_impl;

ifxc_status ifxc_ml25_eval(
    const ifxc_func_type *func,
    const ifxc_input *input,
    size_t nentries,
    const ifxc_deriv_entry *entries);

ifxc_status ifxc_ml25_eval_first_derivatives_contracted(
    const ifxc_handle_impl *impl,
    const ifxc_input *input,
    const double *coeffs,
    double *d_rho,
    double *d_sigma,
    double *d_tau);

ifxc_status ifxc_ml25_eval_integral_order0(
    const ifxc_handle_impl *impl,
    const ifxc_input *input,
    double *features);

ifxc_status ifxc_derivative_component_count(
    const ifxc_dimensions_t *dims,
    const ifxc_deriv_entry *entry,
    size_t *ncomponents);

#endif /* IFXC_INTERNAL_H */
