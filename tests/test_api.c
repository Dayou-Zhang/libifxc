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

static const char *const expected_keys[IFXC_ML25_NFEATURES] = {
#define IFXC_ML25_FEATURE(index, feature_name, feature_key, feature_group, feature_kind) feature_key,
#include "../src/features/ifxc_ml25_features.def"
#undef IFXC_ML25_FEATURE
};

static const char *const expected_groups[IFXC_ML25_NFEATURES] = {
#define IFXC_ML25_FEATURE(index, feature_name, feature_key, feature_group, feature_kind) feature_group,
#include "../src/features/ifxc_ml25_features.def"
#undef IFXC_ML25_FEATURE
};

static const ifxc_feature_kind expected_kinds[IFXC_ML25_NFEATURES] = {
#define IFXC_ML25_FEATURE(index, feature_name, feature_key, feature_group, feature_kind) feature_kind,
#include "../src/features/ifxc_ml25_features.def"
#undef IFXC_ML25_FEATURE
};

static void
check_close(double actual, double expected)
{
  const double abs_tol = 1e-12;
  const double rel_tol = 1e-10;
  double diff = fabs(actual - expected);
  double scale = fabs(expected);

  if(diff > abs_tol && diff > rel_tol * (scale > 1.0 ? scale : 1.0)){
    fprintf(stderr, "mismatch: got %.17g expected %.17g\n", actual, expected);
    assert(diff <= abs_tol || diff <= rel_tol * (scale > 1.0 ? scale : 1.0));
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
  assert(set_info->max_deriv_order == IFXC_MAX_DERIV_ORDER);

  for(i = 0; i < IFXC_ML25_NFEATURES; ++i){
    check_status(ifxc_feature_info(IFXC_FEATURE_SET_ML25, i, &feature_info));
    assert(feature_info != NULL);
    assert(feature_info->index == (int)i);
    assert(feature_info->paper_id == (int)i + 1);
    assert(strcmp(feature_info->key, expected_keys[i]) == 0);
    assert(strcmp(feature_info->group, expected_groups[i]) == 0);
    assert(feature_info->kind == expected_kinds[i]);
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
  assert(max_order == IFXC_MAX_DERIV_ORDER);
  check_status(ifxc_func_dimensions(&func, &dims));
  assert(dims.rho == 2);
  assert(dims.sigma == 3);
  assert(dims.lapl == 0);
  assert(dims.tau == 2);
  ifxc_end(&func);
}

static void
check_integral_matches_local(
    size_t nfeatures,
    size_t npoints,
    const double *local_out,
    const double *integral_out,
    const double *weights)
{
  size_t f;
  size_t p;

  for(f = 0; f < nfeatures; ++f){
    double sum = 0.0;
    for(p = 0; p < npoints; ++p){
      sum += weights[p] * local_out[f * npoints + p];
    }
    check_close(integral_out[f], sum);
  }
}

static void
check_unpolarized_eval(void)
{
  ifxc_func_type func;
  double rho[2] = {0.3, 0.4};
  double sigma[2] = {0.05, 0.02};
  double tau[2] = {0.1, 0.11};
  double weights[2] = {1.0, 0.5};
  double local_out[IFXC_ML25_NFEATURES * 2];
  double integral_out[IFXC_ML25_NFEATURES];
  double rho_deriv_out[IFXC_ML25_NFEATURES * 2];
  double rho_sigma_deriv_out[IFXC_ML25_NFEATURES * 2];
  ifxc_input input = {
    .npoints = 2,
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
  ifxc_deriv_entry rho_deriv_entry = {
    .target = IFXC_TARGET_LOCAL,
    .order = 1,
    .vars = { IFXC_VAR_RHO },
    .out = rho_deriv_out
  };
  ifxc_deriv_entry rho_sigma_deriv_entry = {
    .target = IFXC_TARGET_LOCAL,
    .order = 2,
    .vars = { IFXC_VAR_RHO, IFXC_VAR_SIGMA },
    .out = rho_sigma_deriv_out
  };
  size_t n_double = 0;
  ifxc_status status;

  check_status(ifxc_init(&func, IFXC_FEATURE_SET_ML25, IFXC_UNPOLARIZED));

  status = ifxc_output_size(&func, &input, &local_entry, &n_double);
  check_status(status);
  assert(n_double == (size_t)IFXC_ML25_NFEATURES * 2);
  status = ifxc_output_size(&func, &input, &integral_entry, &n_double);
  check_status(status);
  assert(n_double == IFXC_ML25_NFEATURES);
  status = ifxc_output_size(&func, &input, &rho_deriv_entry, &n_double);
  check_status(status);
  assert(n_double == (size_t)IFXC_ML25_NFEATURES * 2);
  status = ifxc_output_size(&func, &input, &rho_sigma_deriv_entry, &n_double);
  check_status(status);
  assert(n_double == (size_t)IFXC_ML25_NFEATURES * 2);

  status = ifxc_eval(&func, &input, 4, (const ifxc_deriv_entry[]){local_entry, integral_entry, rho_deriv_entry, rho_sigma_deriv_entry});
  check_status(status);

  check_close(local_out[0], -1.1433666468899044);
  check_close(local_out[1], 0.0);
  check_close(local_out[2], -0.07545507559977532);
  check_close(local_out[6], 0.0006086832927199516);

  check_close(integral_out[0], -1.1433666468899044);
  check_close(integral_out[1], -0.07545507559977532);
  check_close(integral_out[2], -0.06465366715086117);
  check_close(integral_out[3], 0.0006086832927199516);

  check_close(rho_deriv_out[0], -1.5584508655337714);
  check_close(rho_deriv_out[2], -0.08802961322078538);

  check_close(rho_sigma_deriv_out[0], 0.26324786692084);
  check_close(rho_sigma_deriv_out[2], -0.11206861044294927);

  check_integral_matches_local(IFXC_ML25_NFEATURES, 2, local_out, integral_out, weights);

  {
    double bad_weights[2] = {1.0, 0.5};
    ifxc_input bad_input = input;
    bad_input.weights = NULL;
    status = ifxc_eval(&func, &bad_input, 1, &integral_entry);
    assert(status == IFXC_E_MISSING_INPUT);
  }

  {
    ifxc_deriv_entry bad_entry = {
      .target = IFXC_TARGET_LOCAL,
      .order = 2,
      .vars = { IFXC_VAR_SIGMA, IFXC_VAR_RHO },
      .out = rho_sigma_deriv_out
    };
    status = ifxc_eval(&func, &input, 1, &bad_entry);
    assert(status == IFXC_E_INVALID_ARGUMENT);
  }

  {
    ifxc_deriv_entry lapl_entry = {
      .target = IFXC_TARGET_LOCAL,
      .order = 1,
      .vars = { IFXC_VAR_LAPL },
      .out = rho_deriv_out
    };
    status = ifxc_output_size(&func, &input, &lapl_entry, &n_double);
    assert(status == IFXC_E_UNSUPPORTED_VARIABLE);
    status = ifxc_eval(&func, &input, 1, &lapl_entry);
    assert(status == IFXC_E_UNSUPPORTED_VARIABLE);
  }

  ifxc_end(&func);
}

static void
check_polarized_eval(void)
{
  ifxc_func_type func;
  double rho[4] = {0.6, 0.7, 0.4, 0.5};
  double sigma[6] = {0.02, 0.03, 0.04, 0.01, 0.015, 0.02};
  double tau[4] = {0.1, 0.12, 0.2, 0.22};
  double weights[2] = {1.0, 0.5};
  double local_out[IFXC_ML25_NFEATURES * 2];
  double integral_out[IFXC_ML25_NFEATURES];
  double rho_deriv_out[IFXC_ML25_NFEATURES * 4];
  double rho_sigma_deriv_out[IFXC_ML25_NFEATURES * 12];
  ifxc_input input = {
    .npoints = 2,
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
  ifxc_deriv_entry rho_deriv_entry = {
    .target = IFXC_TARGET_LOCAL,
    .order = 1,
    .vars = { IFXC_VAR_RHO },
    .out = rho_deriv_out
  };
  ifxc_deriv_entry rho_sigma_deriv_entry = {
    .target = IFXC_TARGET_LOCAL,
    .order = 2,
    .vars = { IFXC_VAR_RHO, IFXC_VAR_SIGMA },
    .out = rho_sigma_deriv_out
  };
  size_t n_double = 0;
  ifxc_status status;

  check_status(ifxc_init(&func, IFXC_FEATURE_SET_ML25, IFXC_POLARIZED));

  status = ifxc_output_size(&func, &input, &local_entry, &n_double);
  check_status(status);
  assert(n_double == (size_t)IFXC_ML25_NFEATURES * 2);
  status = ifxc_output_size(&func, &input, &integral_entry, &n_double);
  check_status(status);
  assert(n_double == IFXC_ML25_NFEATURES);
  status = ifxc_output_size(&func, &input, &rho_deriv_entry, &n_double);
  check_status(status);
  assert(n_double == (size_t)IFXC_ML25_NFEATURES * 4);
  status = ifxc_output_size(&func, &input, &rho_sigma_deriv_entry, &n_double);
  check_status(status);
  assert(n_double == (size_t)IFXC_ML25_NFEATURES * 12);

  status = ifxc_eval(&func, &input, 4, (const ifxc_deriv_entry[]){local_entry, integral_entry, rho_deriv_entry, rho_sigma_deriv_entry});
  check_status(status);

  check_close(local_out[0], -1.73746382183406);
  check_close(local_out[1], 0.0);
  check_close(local_out[2], -0.08200686576627336);
  check_close(local_out[6], 4.590975734960055e-05);

  check_close(integral_out[0], -1.73746382183406);
  check_close(integral_out[1], -0.08200686576627336);
  check_close(integral_out[2], -0.06999673832468495);
  check_close(integral_out[3], 4.590975734960055e-05);

  check_close(rho_deriv_out[0], -2.4592601851436213);
  check_close(rho_deriv_out[2], -2.1668732002232165);

  check_close(rho_sigma_deriv_out[0], 0.01900933973818406);
  check_close(rho_sigma_deriv_out[4], -3.469446951953614e-18);

  check_integral_matches_local(IFXC_ML25_NFEATURES, 2, local_out, integral_out, weights);

  ifxc_end(&func);
}

int
main(void)
{
  check_metadata();
  check_dimensions();
  check_handle();
  check_unpolarized_eval();
  check_polarized_eval();
  return 0;
}
