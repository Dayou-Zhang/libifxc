#include "ifxc.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static void
check_status(ifxc_status status)
{
  if(status != IFXC_OK){
    fprintf(stderr, "unexpected status: %s (%d)\n", ifxc_strerror(status), (int)status);
    assert(status == IFXC_OK);
  }
}

static void
check_metadata(void)
{
  const ifxc_feature_set_info_t *set_info = NULL;
  const ifxc_feature_info_t *feature_info = NULL;
  size_t i;

  check_status(ifxc_feature_set_info(IFXC_FEATURE_SET_ML25, &set_info));
  assert(set_info != NULL);
  assert(set_info->feature_set == IFXC_FEATURE_SET_ML25);
  assert(strcmp(set_info->key, "ml25") == 0);
  assert(strcmp(set_info->name, "ML25 integral features") == 0);
  assert(set_info->nfeatures == IFXC_ML25_NFEATURES);
  assert(set_info->max_deriv_order == 0);

  for(i = 0; i < IFXC_ML25_NFEATURES; ++i){
    check_status(ifxc_feature_info(IFXC_FEATURE_SET_ML25, i, &feature_info));
    assert(feature_info != NULL);
    assert(feature_info->index == (int)i);
    assert(feature_info->paper_id == (int)i + 1);
    assert(feature_info->key != NULL);
    assert(feature_info->group != NULL);
  }
}

static void
check_dimensions(void)
{
  ifxc_dimensions_t dims;

  check_status(ifxc_dimensions(IFXC_FEATURE_SET_ML25, IFXC_UNPOLARIZED, &dims));
  assert(dims.rho == 1);
  assert(dims.sigma == 1);
  assert(dims.lapl == 0);
  assert(dims.tau == 1);

  check_status(ifxc_dimensions(IFXC_FEATURE_SET_ML25, IFXC_POLARIZED, &dims));
  assert(dims.rho == 2);
  assert(dims.sigma == 3);
  assert(dims.lapl == 0);
  assert(dims.tau == 2);
}

static void
check_handle(void)
{
  ifxc_func_type func;
  size_t nfeatures = 0;
  unsigned int max_order = 999;
  ifxc_dimensions_t dims;

  assert(ifxc_version_major() == IFXC_VERSION_MAJOR);
  assert(ifxc_version_minor() == IFXC_VERSION_MINOR);
  assert(ifxc_version_patch() == IFXC_VERSION_PATCH);
  assert(strcmp(ifxc_version_string(), "0.1.0") == 0);

  check_status(ifxc_init(&func, IFXC_FEATURE_SET_ML25, IFXC_POLARIZED));
  check_status(ifxc_nfeatures(&func, &nfeatures));
  assert(nfeatures == IFXC_ML25_NFEATURES);
  check_status(ifxc_max_deriv_order(&func, &max_order));
  assert(max_order == 0);
  check_status(ifxc_func_dimensions(&func, &dims));
  assert(dims.rho == 2);
  assert(dims.sigma == 3);
  assert(dims.lapl == 0);
  assert(dims.tau == 2);
  ifxc_end(&func);
}

static void
check_eval_and_output_size(void)
{
  ifxc_func_type func;
  double rho[2] = {1.0, 2.0};
  double sigma[3] = {0.1, 0.2, 0.3};
  double tau[2] = {0.4, 0.5};
  double weights[1] = {1.0};
  double local_out[IFXC_ML25_NFEATURES];
  double integral_out[IFXC_ML25_NFEATURES];
  ifxc_input input = {
    .npoints = 1,
    .rho = rho,
    .sigma = sigma,
    .lapl = NULL,
    .tau = tau,
    .weights = weights
  };
  ifxc_deriv_entry local_entry = {
    .target = IFXC_TARGET_LOCAL,
    .order = 0,
    .out = local_out
  };
  ifxc_deriv_entry integral_entry = {
    .target = IFXC_TARGET_INTEGRAL,
    .order = 0,
    .out = integral_out
  };
  ifxc_deriv_entry sigma_deriv_entry = {
    .target = IFXC_TARGET_LOCAL,
    .order = 1,
    .vars = { IFXC_VAR_SIGMA },
    .out = local_out
  };
  size_t n_double = 0;

  check_status(ifxc_init(&func, IFXC_FEATURE_SET_ML25, IFXC_POLARIZED));
  check_status(ifxc_output_size(&func, &input, &local_entry, &n_double));
  assert(n_double == IFXC_ML25_NFEATURES);
  check_status(ifxc_output_size(&func, &input, &integral_entry, &n_double));
  assert(n_double == IFXC_ML25_NFEATURES);
  check_status(ifxc_output_size(&func, &input, &sigma_deriv_entry, &n_double));
  assert(n_double == (size_t)IFXC_ML25_NFEATURES * 3);
  check_status(ifxc_eval(&func, &input, 2, (const ifxc_deriv_entry[]){local_entry, integral_entry}));
  ifxc_end(&func);
}

int
main(void)
{
  check_metadata();
  check_dimensions();
  check_handle();
  check_eval_and_output_size();
  return 0;
}
