#include "ifxc_internal.h"

static const ifxc_feature_set_info_t IFXC_ML25_FEATURE_SET_INFO = {
  .feature_set = IFXC_FEATURE_SET_ML25,
  .key = "ml25",
  .name = "ML25 integral features",
  .nfeatures = IFXC_ML25_NFEATURES,
  .max_deriv_order = IFXC_MAX_DERIV_ORDER
};

#define IFXC_ML25_FEATURE(feature_index, feature_name, feature_key, feature_group, feature_kind) \
  { .index = (feature_index), .paper_id = (feature_index) + 1, .key = (feature_key), .group = (feature_group), .kind = (feature_kind) },

static const ifxc_feature_info_t IFXC_ML25_FEATURES[IFXC_ML25_NFEATURES] = {
#include "features/ifxc_ml25_features.def"
};

#undef IFXC_ML25_FEATURE

ifxc_status
ifxc_feature_set_info(int feature_set, const ifxc_feature_set_info_t **info)
{
  if(info == NULL){
    return IFXC_E_INVALID_ARGUMENT;
  }

  switch(feature_set){
  case IFXC_FEATURE_SET_ML25:
    *info = &IFXC_ML25_FEATURE_SET_INFO;
    return IFXC_OK;
  default:
    *info = NULL;
    return IFXC_E_UNKNOWN_FEATURE_SET;
  }
}

ifxc_status
ifxc_feature_info(int feature_set, size_t index, const ifxc_feature_info_t **info)
{
  if(info == NULL){
    return IFXC_E_INVALID_ARGUMENT;
  }

  switch(feature_set){
  case IFXC_FEATURE_SET_ML25:
    if(index >= IFXC_ML25_NFEATURES){
      *info = NULL;
      return IFXC_E_INVALID_ARGUMENT;
    }
    *info = &IFXC_ML25_FEATURES[index];
    return IFXC_OK;
  default:
    *info = NULL;
    return IFXC_E_UNKNOWN_FEATURE_SET;
  }
}

ifxc_status
ifxc_dimensions(int feature_set, ifxc_nspin nspin, ifxc_dimensions_t *dims)
{
  if(dims == NULL){
    return IFXC_E_INVALID_ARGUMENT;
  }
  if(feature_set != IFXC_FEATURE_SET_ML25){
    return IFXC_E_UNKNOWN_FEATURE_SET;
  }

  switch(nspin){
  case IFXC_UNPOLARIZED:
    dims->rho = 1;
    dims->sigma = 1;
    dims->lapl = 0;
    dims->tau = 1;
    return IFXC_OK;
  case IFXC_POLARIZED:
    dims->rho = 2;
    dims->sigma = 3;
    dims->lapl = 0;
    dims->tau = 2;
    return IFXC_OK;
  default:
    return IFXC_E_INVALID_NSPIN;
  }
}
