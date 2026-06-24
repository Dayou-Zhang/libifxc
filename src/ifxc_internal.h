#ifndef IFXC_INTERNAL_H
#define IFXC_INTERNAL_H

#include "ifxc.h"

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

#endif /* IFXC_INTERNAL_H */
