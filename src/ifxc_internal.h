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

#endif /* IFXC_INTERNAL_H */
