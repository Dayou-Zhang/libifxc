#include "ifxc.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

enum { TEST_ML25_MAX_ORDER = 3 };

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

static const int expected_public_indices[IFXC_ML25_NFEATURES] = {
#define IFXC_ML25_FEATURE(index, feature_name, feature_key, feature_group, feature_kind) IFXC_ML25_##feature_name,
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

static size_t
local_index(size_t nfeatures, size_t point, size_t feature)
{
  return point * nfeatures + feature;
}

static size_t
deriv_index(size_t nfeatures, size_t ncomp, size_t point, size_t feature, size_t comp)
{
  return (point * nfeatures + feature) * ncomp + comp;
}

static size_t
binomial(size_t n, size_t k)
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

static size_t
variable_dimension(const ifxc_dimensions_t *dims, ifxc_variable var)
{
  switch(var){
  case IFXC_VAR_RHO: return dims->rho;
  case IFXC_VAR_SIGMA: return dims->sigma;
  case IFXC_VAR_LAPL: return dims->lapl;
  case IFXC_VAR_TAU: return dims->tau;
  default: return 0;
  }
}

static size_t
expected_component_count(const ifxc_dimensions_t *dims,
                         unsigned int order,
                         const ifxc_variable *vars)
{
  unsigned int counts[4] = {0, 0, 0, 0};
  unsigned int i;
  size_t total = 1;

  if(order == 0){
    return 1;
  }

  for(i = 0; i < order; ++i){
    counts[vars[i]] += 1;
  }
  for(i = 0; i <= (unsigned int)IFXC_VAR_TAU; ++i){
    if(counts[i] > 0){
      size_t dim = variable_dimension(dims, (ifxc_variable)i);
      total *= binomial(dim + counts[i] - 1, counts[i]);
    }
  }

  return total;
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
  assert(set_info->max_deriv_order == TEST_ML25_MAX_ORDER);

  for(i = 0; i < IFXC_ML25_NFEATURES; ++i){
    check_status(ifxc_feature_info(IFXC_FEATURE_SET_ML25, i, &feature_info));
    assert(feature_info != NULL);
    assert(feature_info->index == (int)i);
    assert(expected_public_indices[i] == (int)i);
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
  assert(max_order == TEST_ML25_MAX_ORDER);
  check_status(ifxc_func_dimensions(&func, &dims));
  assert(dims.rho == 2);
  assert(dims.sigma == 3);
  assert(dims.lapl == 0);
  assert(dims.tau == 2);
  ifxc_end(&func);
}

static void
check_error_paths(void)
{
  ifxc_func_type func;
  ifxc_func_type uninitialized = {0};
  ifxc_dimensions_t dims;
  const ifxc_feature_set_info_t *set_info = NULL;
  const ifxc_feature_info_t *feature_info = NULL;
  double rho[1] = {0.3};
  double sigma[1] = {0.05};
  double tau[1] = {0.1};
  double out[IFXC_ML25_NFEATURES];
  ifxc_input input = {
    .npoints = 1,
    .rho = rho,
    .sigma = sigma,
    .lapl = NULL,
    .tau = tau,
    .weights = NULL
  };
  ifxc_deriv_entry entry = {
    .target = IFXC_TARGET_LOCAL,
    .order = 0,
    .out = out
  };
  size_t nfeatures = 0;
  unsigned int max_order = 0;
  size_t n_double = 0;

  assert(ifxc_init(NULL, IFXC_FEATURE_SET_ML25, IFXC_UNPOLARIZED) == IFXC_E_INVALID_ARGUMENT);
  assert(ifxc_init(&func, 999, IFXC_UNPOLARIZED) == IFXC_E_UNKNOWN_FEATURE_SET);
  assert(ifxc_init(&func, IFXC_FEATURE_SET_ML25, (ifxc_nspin)99) == IFXC_E_INVALID_NSPIN);
  assert(ifxc_dimensions(999, IFXC_UNPOLARIZED, &dims) == IFXC_E_UNKNOWN_FEATURE_SET);
  assert(ifxc_dimensions(IFXC_FEATURE_SET_ML25, (ifxc_nspin)99, &dims) == IFXC_E_INVALID_NSPIN);
  assert(ifxc_dimensions(IFXC_FEATURE_SET_ML25, IFXC_UNPOLARIZED, NULL) == IFXC_E_INVALID_ARGUMENT);
  assert(ifxc_feature_set_info(999, &set_info) == IFXC_E_UNKNOWN_FEATURE_SET);
  assert(ifxc_feature_set_info(IFXC_FEATURE_SET_ML25, NULL) == IFXC_E_INVALID_ARGUMENT);
  assert(ifxc_feature_info(IFXC_FEATURE_SET_ML25, IFXC_ML25_NFEATURES, &feature_info) == IFXC_E_INVALID_ARGUMENT);
  assert(ifxc_feature_info(IFXC_FEATURE_SET_ML25, 0, NULL) == IFXC_E_INVALID_ARGUMENT);

  assert(ifxc_nfeatures(&uninitialized, &nfeatures) == IFXC_E_NOT_INITIALIZED);
  assert(ifxc_max_deriv_order(&uninitialized, &max_order) == IFXC_E_NOT_INITIALIZED);
  assert(ifxc_func_dimensions(&uninitialized, &dims) == IFXC_E_NOT_INITIALIZED);
  assert(ifxc_output_size(&uninitialized, &input, &entry, &n_double) == IFXC_E_NOT_INITIALIZED);
  assert(ifxc_eval(&uninitialized, &input, 1, &entry) == IFXC_E_NOT_INITIALIZED);

  check_status(ifxc_init(&func, IFXC_FEATURE_SET_ML25, IFXC_UNPOLARIZED));
  assert(ifxc_nfeatures(&func, NULL) == IFXC_E_INVALID_ARGUMENT);
  assert(ifxc_max_deriv_order(&func, NULL) == IFXC_E_INVALID_ARGUMENT);
  assert(ifxc_func_dimensions(&func, NULL) == IFXC_E_INVALID_ARGUMENT);
  assert(ifxc_output_size(&func, &input, &entry, NULL) == IFXC_E_INVALID_ARGUMENT);
  assert(ifxc_eval(&func, &input, 0, &entry) == IFXC_E_INVALID_ARGUMENT);
  assert(ifxc_eval(&func, &input, 1, NULL) == IFXC_E_INVALID_ARGUMENT);

  {
    ifxc_deriv_entry sizing_entry = entry;
    sizing_entry.out = NULL;
    check_status(ifxc_output_size(&func, &input, &sizing_entry, &n_double));
    assert(n_double == IFXC_ML25_NFEATURES * input.npoints);
    assert(ifxc_eval(&func, &input, 1, &sizing_entry) == IFXC_E_INVALID_ARGUMENT);
  }

  {
    ifxc_input missing_rho = input;
    ifxc_input missing_sigma = input;
    ifxc_input missing_tau = input;
    missing_rho.rho = NULL;
    missing_sigma.sigma = NULL;
    missing_tau.tau = NULL;
    assert(ifxc_eval(&func, &missing_rho, 1, &entry) == IFXC_E_MISSING_INPUT);
    assert(ifxc_eval(&func, &missing_sigma, 1, &entry) == IFXC_E_MISSING_INPUT);
    assert(ifxc_eval(&func, &missing_tau, 1, &entry) == IFXC_E_MISSING_INPUT);
  }

  ifxc_end(&func);
}

static void
check_output_size_matrix_for_spin(ifxc_nspin nspin)
{
  static const ifxc_variable vars_order1[][1] = {
    { IFXC_VAR_RHO },
    { IFXC_VAR_SIGMA },
    { IFXC_VAR_TAU }
  };
  static const ifxc_variable vars_order2[][2] = {
    { IFXC_VAR_RHO, IFXC_VAR_RHO },
    { IFXC_VAR_RHO, IFXC_VAR_SIGMA },
    { IFXC_VAR_RHO, IFXC_VAR_TAU },
    { IFXC_VAR_SIGMA, IFXC_VAR_SIGMA },
    { IFXC_VAR_SIGMA, IFXC_VAR_TAU },
    { IFXC_VAR_TAU, IFXC_VAR_TAU }
  };
  static const ifxc_variable vars_order3[][3] = {
    { IFXC_VAR_RHO, IFXC_VAR_RHO, IFXC_VAR_RHO },
    { IFXC_VAR_RHO, IFXC_VAR_RHO, IFXC_VAR_SIGMA },
    { IFXC_VAR_RHO, IFXC_VAR_RHO, IFXC_VAR_TAU },
    { IFXC_VAR_RHO, IFXC_VAR_SIGMA, IFXC_VAR_SIGMA },
    { IFXC_VAR_RHO, IFXC_VAR_SIGMA, IFXC_VAR_TAU },
    { IFXC_VAR_RHO, IFXC_VAR_TAU, IFXC_VAR_TAU },
    { IFXC_VAR_SIGMA, IFXC_VAR_SIGMA, IFXC_VAR_SIGMA },
    { IFXC_VAR_SIGMA, IFXC_VAR_SIGMA, IFXC_VAR_TAU },
    { IFXC_VAR_SIGMA, IFXC_VAR_TAU, IFXC_VAR_TAU },
    { IFXC_VAR_TAU, IFXC_VAR_TAU, IFXC_VAR_TAU }
  };
  ifxc_func_type func;
  ifxc_dimensions_t dims;
  double rho[8] = {0};
  double sigma[12] = {0};
  double tau[8] = {0};
  ifxc_input input = {
    .npoints = 4,
    .rho = rho,
    .sigma = sigma,
    .lapl = NULL,
    .tau = tau,
    .weights = NULL
  };
  size_t n_double = 0;
  size_t i;
  ifxc_deriv_entry entry = {
    .target = IFXC_TARGET_LOCAL,
    .order = 0,
    .vars = NULL,
    .out = NULL
  };

  check_status(ifxc_init(&func, IFXC_FEATURE_SET_ML25, nspin));
  check_status(ifxc_func_dimensions(&func, &dims));

  check_status(ifxc_output_size(&func, &input, &entry, &n_double));
  assert(n_double == (size_t)IFXC_ML25_NFEATURES * input.npoints);
  entry.target = IFXC_TARGET_INTEGRAL;
  check_status(ifxc_output_size(&func, &input, &entry, &n_double));
  assert(n_double == IFXC_ML25_NFEATURES);
  entry.target = IFXC_TARGET_LOCAL;

  for(i = 0; i < sizeof(vars_order1) / sizeof(vars_order1[0]); ++i){
    entry.order = 1;
    entry.vars = vars_order1[i];
    check_status(ifxc_output_size(&func, &input, &entry, &n_double));
    assert(n_double == (size_t)IFXC_ML25_NFEATURES * input.npoints *
                       expected_component_count(&dims, entry.order, entry.vars));
  }
  for(i = 0; i < sizeof(vars_order2) / sizeof(vars_order2[0]); ++i){
    entry.order = 2;
    entry.vars = vars_order2[i];
    check_status(ifxc_output_size(&func, &input, &entry, &n_double));
    assert(n_double == (size_t)IFXC_ML25_NFEATURES * input.npoints *
                       expected_component_count(&dims, entry.order, entry.vars));
  }
  for(i = 0; i < sizeof(vars_order3) / sizeof(vars_order3[0]); ++i){
    entry.order = 3;
    entry.vars = vars_order3[i];
    check_status(ifxc_output_size(&func, &input, &entry, &n_double));
    assert(n_double == (size_t)IFXC_ML25_NFEATURES * input.npoints *
                       expected_component_count(&dims, entry.order, entry.vars));
  }

  ifxc_end(&func);
}

static void
check_output_size_matrix(void)
{
  check_output_size_matrix_for_spin(IFXC_UNPOLARIZED);
  check_output_size_matrix_for_spin(IFXC_POLARIZED);
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
      sum += weights[p] * local_out[local_index(nfeatures, p, f)];
    }
    check_close(integral_out[f], sum);
  }
}

static void
check_weighted_deriv_matches_local(
    size_t nfeatures,
    size_t npoints,
    size_t ncomp,
    const double *local_out,
    const double *weighted_out,
    const double *weights)
{
  size_t f;
  size_t p;
  size_t c;

  for(p = 0; p < npoints; ++p){
    for(f = 0; f < nfeatures; ++f){
      for(c = 0; c < ncomp; ++c){
        size_t idx = deriv_index(nfeatures, ncomp, p, f, c);
        check_close(weighted_out[idx], weights[p] * local_out[idx]);
      }
    }
  }
}

static void
check_ml25_sample_feature_values(
    size_t nfeatures,
    size_t point,
    const double *local_out,
    const double *expected_values)
{
  static const size_t features[] = {
    IFXC_ML25_LAK_X,
    IFXC_ML25_LAK_C,
    IFXC_ML25_LYP_T1,
    IFXC_ML25_MN15_C01,
    IFXC_ML25_MN15_C20,
    IFXC_ML25_MN15_A01,
    IFXC_ML25_MN15_B01,
    IFXC_ML25_MN15_B09
  };
  size_t i;

  for(i = 0; i < sizeof(features) / sizeof(features[0]); ++i){
    check_close(local_out[local_index(nfeatures, point, features[i])],
                expected_values[i]);
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
  double rho_integral_deriv_out[IFXC_ML25_NFEATURES * 2];
  double rho_sigma_deriv_out[IFXC_ML25_NFEATURES * 2];
  double rho_rho_tau_deriv_out[IFXC_ML25_NFEATURES * 2];
  const ifxc_variable rho_vars[] = { IFXC_VAR_RHO };
  const ifxc_variable rho_sigma_vars[] = { IFXC_VAR_RHO, IFXC_VAR_SIGMA };
  const ifxc_variable rho_rho_tau_vars[] = { IFXC_VAR_RHO, IFXC_VAR_RHO, IFXC_VAR_TAU };
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
    .vars = rho_vars,
    .out = rho_deriv_out
  };
  ifxc_deriv_entry rho_integral_deriv_entry = {
    .target = IFXC_TARGET_INTEGRAL,
    .order = 1,
    .vars = rho_vars,
    .out = rho_integral_deriv_out
  };
  ifxc_deriv_entry rho_sigma_deriv_entry = {
    .target = IFXC_TARGET_LOCAL,
    .order = 2,
    .vars = rho_sigma_vars,
    .out = rho_sigma_deriv_out
  };
  ifxc_deriv_entry rho_rho_tau_deriv_entry = {
    .target = IFXC_TARGET_LOCAL,
    .order = 3,
    .vars = rho_rho_tau_vars,
    .out = rho_rho_tau_deriv_out
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
  status = ifxc_output_size(&func, &input, &rho_integral_deriv_entry, &n_double);
  check_status(status);
  assert(n_double == (size_t)IFXC_ML25_NFEATURES * 2);
  status = ifxc_output_size(&func, &input, &rho_sigma_deriv_entry, &n_double);
  check_status(status);
  assert(n_double == (size_t)IFXC_ML25_NFEATURES * 2);
  status = ifxc_output_size(&func, &input, &rho_rho_tau_deriv_entry, &n_double);
  check_status(status);
  assert(n_double == (size_t)IFXC_ML25_NFEATURES * 2);
  status = ifxc_eval(&func, &input, 6, (const ifxc_deriv_entry[]){
    local_entry,
    integral_entry,
    rho_deriv_entry,
    rho_integral_deriv_entry,
    rho_sigma_deriv_entry,
    rho_rho_tau_deriv_entry
  });
  check_status(status);

  check_close(local_out[local_index(IFXC_ML25_NFEATURES, 0, IFXC_ML25_LAK_X)],
              -0.3430099940669713);
  check_close(local_out[local_index(IFXC_ML25_NFEATURES, 0, IFXC_ML25_LAK_C)],
              -0.022636522679932596);
  check_close(local_out[local_index(IFXC_ML25_NFEATURES, 0, IFXC_ML25_LYP_T1)],
              -0.01939610014525835);
  check_close(local_out[local_index(IFXC_ML25_NFEATURES, 0, IFXC_ML25_LYP_T2)],
              0.00018260498781598548);
  {
    static const double expected[] = {
      -0.3430099940669713,
      -0.022636522679932596,
      -0.01939610014525835,
      -0.29664934427289896,
      -0.09959623703114227,
      -0.03692412583283123,
      0.001997182896612906,
      2.8729640128347612e-05
    };
    check_ml25_sample_feature_values(IFXC_ML25_NFEATURES, 0, local_out, expected);
  }
  assert(fabs(local_out[local_index(IFXC_ML25_NFEATURES, 1, IFXC_ML25_LAK_X)]) > 1e-12);
  assert(fabs(local_out[local_index(IFXC_ML25_NFEATURES, 1, IFXC_ML25_LAK_C)]) > 1e-12);

  check_close(rho_deriv_out[deriv_index(IFXC_ML25_NFEATURES, 1, 0, IFXC_ML25_LAK_X, 0)],
              -1.5584508655337714);
  check_close(rho_deriv_out[deriv_index(IFXC_ML25_NFEATURES, 1, 0, IFXC_ML25_LAK_C, 0)],
              -0.08802961322078538);

  check_close(rho_sigma_deriv_out[deriv_index(IFXC_ML25_NFEATURES, 1, 0, IFXC_ML25_LAK_X, 0)],
              0.26324786692084);
  check_close(rho_sigma_deriv_out[deriv_index(IFXC_ML25_NFEATURES, 1, 0, IFXC_ML25_LAK_C, 0)],
              -0.11206861044294927);

  assert(isfinite(rho_rho_tau_deriv_out[0]));

  check_integral_matches_local(IFXC_ML25_NFEATURES, 2, local_out, integral_out, weights);
  check_weighted_deriv_matches_local(IFXC_ML25_NFEATURES, 2, 1,
                                     rho_deriv_out, rho_integral_deriv_out, weights);

  {
    ifxc_input bad_input = input;
    bad_input.weights = NULL;
    status = ifxc_eval(&func, &bad_input, 1, &integral_entry);
    assert(status == IFXC_E_MISSING_INPUT);
  }

  {
    const ifxc_variable bad_vars[] = { IFXC_VAR_SIGMA, IFXC_VAR_RHO };
    ifxc_deriv_entry bad_entry = {
      .target = IFXC_TARGET_LOCAL,
      .order = 2,
      .vars = bad_vars,
      .out = rho_sigma_deriv_out
    };
    status = ifxc_eval(&func, &input, 1, &bad_entry);
    assert(status == IFXC_E_INVALID_ARGUMENT);
  }

  {
    const ifxc_variable lapl_vars[] = { IFXC_VAR_LAPL };
    ifxc_deriv_entry lapl_entry = {
      .target = IFXC_TARGET_LOCAL,
      .order = 1,
      .vars = lapl_vars,
      .out = rho_deriv_out
    };
    status = ifxc_output_size(&func, &input, &lapl_entry, &n_double);
    assert(status == IFXC_E_UNSUPPORTED_VARIABLE);
    status = ifxc_eval(&func, &input, 1, &lapl_entry);
    assert(status == IFXC_E_UNSUPPORTED_VARIABLE);
  }

  {
    const ifxc_variable order4_vars[] = {
      IFXC_VAR_RHO, IFXC_VAR_RHO, IFXC_VAR_TAU, IFXC_VAR_TAU
    };
    ifxc_deriv_entry unsupported_entry = {
      .target = IFXC_TARGET_LOCAL,
      .order = 4,
      .vars = order4_vars,
      .out = rho_deriv_out
    };
    status = ifxc_output_size(&func, &input, &unsupported_entry, &n_double);
    assert(status == IFXC_E_UNSUPPORTED_DERIVATIVE);
    status = ifxc_eval(&func, &input, 1, &unsupported_entry);
    assert(status == IFXC_E_UNSUPPORTED_DERIVATIVE);
  }

  {
    ifxc_deriv_entry missing_vars_entry = {
      .target = IFXC_TARGET_LOCAL,
      .order = 1,
      .vars = NULL,
      .out = rho_deriv_out
    };
    status = ifxc_eval(&func, &input, 1, &missing_vars_entry);
    assert(status == IFXC_E_INVALID_ARGUMENT);
  }

  {
    ifxc_deriv_entry bad_target_entry = rho_deriv_entry;
    bad_target_entry.target = (ifxc_target)99;
    status = ifxc_output_size(&func, &input, &bad_target_entry, &n_double);
    assert(status == IFXC_E_INVALID_ARGUMENT);
    status = ifxc_eval(&func, &input, 1, &bad_target_entry);
    assert(status == IFXC_E_INVALID_ARGUMENT);
  }

  ifxc_end(&func);
}

static void
check_unpolarized_order0_only_eval(void)
{
  ifxc_func_type func;
  double rho[2] = {0.3, 0.4};
  double sigma[2] = {0.05, 0.02};
  double tau[2] = {0.1, 0.11};
  double weights[2] = {1.0, 0.5};
  double local_out[IFXC_ML25_NFEATURES * 2];
  double integral_out[IFXC_ML25_NFEATURES];
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

  check_status(ifxc_init(&func, IFXC_FEATURE_SET_ML25, IFXC_UNPOLARIZED));
  check_status(ifxc_eval(&func, &input, 2, (const ifxc_deriv_entry[]){
    local_entry,
    integral_entry
  }));

  check_close(local_out[local_index(IFXC_ML25_NFEATURES, 0, IFXC_ML25_LAK_X)],
              -0.3430099940669713);
  check_close(local_out[local_index(IFXC_ML25_NFEATURES, 0, IFXC_ML25_LAK_C)],
              -0.022636522679932596);
  assert(fabs(local_out[local_index(IFXC_ML25_NFEATURES, 1, IFXC_ML25_LAK_X)]) > 1e-12);
  check_integral_matches_local(IFXC_ML25_NFEATURES, 2, local_out, integral_out, weights);

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
  double rho_deriv_out[IFXC_ML25_NFEATURES * 2 * 2];
  double rho_integral_deriv_out[IFXC_ML25_NFEATURES * 2 * 2];
  double rho_sigma_deriv_out[IFXC_ML25_NFEATURES * 2 * 6];
  double rho_rho_tau_deriv_out[IFXC_ML25_NFEATURES * 2 * 6];
  const ifxc_variable rho_vars[] = { IFXC_VAR_RHO };
  const ifxc_variable rho_sigma_vars[] = { IFXC_VAR_RHO, IFXC_VAR_SIGMA };
  const ifxc_variable rho_rho_tau_vars[] = { IFXC_VAR_RHO, IFXC_VAR_RHO, IFXC_VAR_TAU };
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
    .vars = rho_vars,
    .out = rho_deriv_out
  };
  ifxc_deriv_entry rho_integral_deriv_entry = {
    .target = IFXC_TARGET_INTEGRAL,
    .order = 1,
    .vars = rho_vars,
    .out = rho_integral_deriv_out
  };
  ifxc_deriv_entry rho_sigma_deriv_entry = {
    .target = IFXC_TARGET_LOCAL,
    .order = 2,
    .vars = rho_sigma_vars,
    .out = rho_sigma_deriv_out
  };
  ifxc_deriv_entry rho_rho_tau_deriv_entry = {
    .target = IFXC_TARGET_LOCAL,
    .order = 3,
    .vars = rho_rho_tau_vars,
    .out = rho_rho_tau_deriv_out
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
  assert(n_double == (size_t)IFXC_ML25_NFEATURES * 2 * 2);
  status = ifxc_output_size(&func, &input, &rho_integral_deriv_entry, &n_double);
  check_status(status);
  assert(n_double == (size_t)IFXC_ML25_NFEATURES * 2 * 2);
  status = ifxc_output_size(&func, &input, &rho_sigma_deriv_entry, &n_double);
  check_status(status);
  assert(n_double == (size_t)IFXC_ML25_NFEATURES * 2 * 6);
  status = ifxc_output_size(&func, &input, &rho_rho_tau_deriv_entry, &n_double);
  check_status(status);
  assert(n_double == (size_t)IFXC_ML25_NFEATURES * 2 * 6);
  status = ifxc_eval(&func, &input, 6, (const ifxc_deriv_entry[]){
    local_entry,
    integral_entry,
    rho_deriv_entry,
    rho_integral_deriv_entry,
    rho_sigma_deriv_entry,
    rho_rho_tau_deriv_entry
  });
  check_status(status);

  check_close(local_out[local_index(IFXC_ML25_NFEATURES, 0, IFXC_ML25_LAK_X)],
              -1.73746382183406);
  check_close(local_out[local_index(IFXC_ML25_NFEATURES, 0, IFXC_ML25_LAK_C)],
              -0.08200686576627336);
  check_close(local_out[local_index(IFXC_ML25_NFEATURES, 0, IFXC_ML25_LYP_T1)],
              -0.06999673832468495);
  check_close(local_out[local_index(IFXC_ML25_NFEATURES, 0, IFXC_ML25_LYP_T2)],
              4.590975734960055e-05);
  assert(fabs(local_out[local_index(IFXC_ML25_NFEATURES, 1, IFXC_ML25_LAK_X)]) > 1e-12);
  assert(fabs(local_out[local_index(IFXC_ML25_NFEATURES, 1, IFXC_ML25_LAK_C)]) > 1e-12);

  check_close(rho_deriv_out[deriv_index(IFXC_ML25_NFEATURES, 2, 0, IFXC_ML25_LAK_X, 0)],
              -2.4592601851436213);
  check_close(rho_deriv_out[deriv_index(IFXC_ML25_NFEATURES, 2, 0, IFXC_ML25_LAK_X, 1)],
              -2.1668732002232165);

  check_close(rho_sigma_deriv_out[deriv_index(IFXC_ML25_NFEATURES, 6, 0, IFXC_ML25_LAK_X, 0)],
              0.01900933973818406);
  check_close(rho_sigma_deriv_out[deriv_index(IFXC_ML25_NFEATURES, 6, 0, IFXC_ML25_LAK_X, 2)],
              -3.469446951953614e-18);

  assert(isfinite(rho_rho_tau_deriv_out[0]));

  check_integral_matches_local(IFXC_ML25_NFEATURES, 2, local_out, integral_out, weights);
  check_weighted_deriv_matches_local(IFXC_ML25_NFEATURES, 2, 2,
                                     rho_deriv_out, rho_integral_deriv_out, weights);

  ifxc_end(&func);
}

static void
check_polarized_ml25_sample_values(void)
{
  ifxc_func_type func;
  double rho[2] = {0.30, 0.20};
  double sigma[3] = {0.06, 0.03, 0.05};
  double tau[2] = {0.10, 0.08};
  double local_out[IFXC_ML25_NFEATURES];
  ifxc_input input = {
    .npoints = 1,
    .rho = rho,
    .sigma = sigma,
    .lapl = NULL,
    .tau = tau,
    .weights = NULL
  };
  ifxc_deriv_entry entry = {
    .target = IFXC_TARGET_LOCAL,
    .order = 0,
    .out = local_out
  };
  static const double expected[] = {
    -0.6882215860668606,
    -0.03785794239322413,
    -0.03279321514349554,
    -0.5914246540028052,
    -0.2447765626418307,
    -0.06468940177494086,
    0.0034137750998256312,
    0.0001353685476186051
  };

  check_status(ifxc_init(&func, IFXC_FEATURE_SET_ML25, IFXC_POLARIZED));
  check_status(ifxc_eval(&func, &input, 1, &entry));
  check_ml25_sample_feature_values(IFXC_ML25_NFEATURES, 0, local_out, expected);
  ifxc_end(&func);
}

static void
check_deterministic_grid_for_spin(ifxc_nspin nspin)
{
  ifxc_func_type func;
  size_t npoints = 4;
  double weights[4] = {0.7, 1.1, 0.4, 0.9};
  double local_out[IFXC_ML25_NFEATURES * 4];
  double integral_out[IFXC_ML25_NFEATURES];
  double rho_unpol[4] = {0.28, 0.36, 0.44, 0.52};
  double sigma_unpol[4] = {0.018, 0.026, 0.034, 0.042};
  double tau_unpol[4] = {0.08, 0.10, 0.13, 0.16};
  double rho_pol[8] = {
    0.33, 0.41, 0.49, 0.57,
    0.24, 0.31, 0.38, 0.45
  };
  double sigma_pol[12] = {
    0.015, 0.019, 0.024, 0.030,
    0.006, 0.008, 0.011, 0.014,
    0.017, 0.023, 0.028, 0.035
  };
  double tau_pol[8] = {
    0.070, 0.090, 0.115, 0.140,
    0.060, 0.080, 0.100, 0.125
  };
  ifxc_input input = {
    .npoints = npoints,
    .rho = (nspin == IFXC_UNPOLARIZED) ? rho_unpol : rho_pol,
    .sigma = (nspin == IFXC_UNPOLARIZED) ? sigma_unpol : sigma_pol,
    .lapl = NULL,
    .tau = (nspin == IFXC_UNPOLARIZED) ? tau_unpol : tau_pol,
    .weights = weights
  };
  ifxc_deriv_entry entries[2] = {
    {
      .target = IFXC_TARGET_LOCAL,
      .order = 0,
      .out = local_out
    },
    {
      .target = IFXC_TARGET_INTEGRAL,
      .order = 0,
      .out = integral_out
    }
  };
  size_t p;

  check_status(ifxc_init(&func, IFXC_FEATURE_SET_ML25, nspin));
  check_status(ifxc_eval(&func, &input, 2, entries));

  for(p = 0; p < npoints; ++p){
    assert(fabs(local_out[local_index(IFXC_ML25_NFEATURES, p, IFXC_ML25_LAK_X)]) > 1e-12);
    assert(fabs(local_out[local_index(IFXC_ML25_NFEATURES, p, IFXC_ML25_LAK_C)]) > 1e-12);
    assert(fabs(local_out[local_index(IFXC_ML25_NFEATURES, p, IFXC_ML25_LYP_T1)]) > 1e-12);
    assert(fabs(local_out[local_index(IFXC_ML25_NFEATURES, p, IFXC_ML25_MN15_C20)]) > 1e-12);
    assert(fabs(local_out[local_index(IFXC_ML25_NFEATURES, p, IFXC_ML25_MN15_A01)]) > 1e-12);
    assert(fabs(local_out[local_index(IFXC_ML25_NFEATURES, p, IFXC_ML25_MN15_B01)]) > 1e-12);
  }

  check_integral_matches_local(IFXC_ML25_NFEATURES, npoints, local_out, integral_out, weights);
  ifxc_end(&func);
}

static void
eval_unpolarized_feature_value(size_t feature,
                               double rho,
                               double sigma,
                               double tau,
                               double *value)
{
  ifxc_func_type func;
  double local_out[IFXC_ML25_NFEATURES];
  ifxc_input input = {
    .npoints = 1,
    .rho = &rho,
    .sigma = &sigma,
    .lapl = NULL,
    .tau = &tau,
    .weights = NULL
  };
  ifxc_deriv_entry entry = {
    .target = IFXC_TARGET_LOCAL,
    .order = 0,
    .out = local_out
  };

  check_status(ifxc_init(&func, IFXC_FEATURE_SET_ML25, IFXC_UNPOLARIZED));
  check_status(ifxc_eval(&func, &input, 1, &entry));
  *value = local_out[feature];
  ifxc_end(&func);
}

static void
eval_unpolarized_feature_derivatives(double rho,
                                     double sigma,
                                     double tau,
                                     double *drho,
                                     double *dsigma,
                                     double *dtau)
{
  ifxc_func_type func;
  double rho_out[IFXC_ML25_NFEATURES];
  double sigma_out[IFXC_ML25_NFEATURES];
  double tau_out[IFXC_ML25_NFEATURES];
  const ifxc_variable rho_vars[] = { IFXC_VAR_RHO };
  const ifxc_variable sigma_vars[] = { IFXC_VAR_SIGMA };
  const ifxc_variable tau_vars[] = { IFXC_VAR_TAU };
  ifxc_input input = {
    .npoints = 1,
    .rho = &rho,
    .sigma = &sigma,
    .lapl = NULL,
    .tau = &tau,
    .weights = NULL
  };
  ifxc_deriv_entry entries[3] = {
    {
      .target = IFXC_TARGET_LOCAL,
      .order = 1,
      .vars = rho_vars,
      .out = rho_out
    },
    {
      .target = IFXC_TARGET_LOCAL,
      .order = 1,
      .vars = sigma_vars,
      .out = sigma_out
    },
    {
      .target = IFXC_TARGET_LOCAL,
      .order = 1,
      .vars = tau_vars,
      .out = tau_out
    }
  };

  check_status(ifxc_init(&func, IFXC_FEATURE_SET_ML25, IFXC_UNPOLARIZED));
  check_status(ifxc_eval(&func, &input, 3, entries));
  memcpy(drho, rho_out, sizeof(rho_out));
  memcpy(dsigma, sigma_out, sizeof(sigma_out));
  memcpy(dtau, tau_out, sizeof(tau_out));
  ifxc_end(&func);
}

static void
check_finite_difference_first_derivatives(void)
{
  static const size_t features[] = {
    IFXC_ML25_LAK_X,
    IFXC_ML25_LAK_C,
    IFXC_ML25_LYP_T1,
    IFXC_ML25_LYP_T4,
    IFXC_ML25_MN15_C01,
    IFXC_ML25_MN15_C20,
    IFXC_ML25_MN15_A01,
    IFXC_ML25_MN15_B01
  };
  const double rho = 0.42;
  const double sigma = 0.031;
  const double tau = 0.12;
  const double h_rho = 1e-5;
  const double h_sigma = 1e-6;
  const double h_tau = 1e-6;
  double drho[IFXC_ML25_NFEATURES];
  double dsigma[IFXC_ML25_NFEATURES];
  double dtau[IFXC_ML25_NFEATURES];
  size_t i;

  eval_unpolarized_feature_derivatives(rho, sigma, tau, drho, dsigma, dtau);

  for(i = 0; i < sizeof(features) / sizeof(features[0]); ++i){
    size_t feature = features[i];
    double plus;
    double minus;
    double fd;

    eval_unpolarized_feature_value(feature, rho + h_rho, sigma, tau, &plus);
    eval_unpolarized_feature_value(feature, rho - h_rho, sigma, tau, &minus);
    fd = (plus - minus) / (2.0 * h_rho);
    if(fabs(drho[feature] - fd) > 1e-8 * (fabs(fd) > 1.0 ? fabs(fd) : 1.0)){
      fprintf(stderr, "finite difference rho feature %zu analytic %.17g fd %.17g plus %.17g minus %.17g\n",
              feature, drho[feature], fd, plus, minus);
    }
    check_close(drho[feature], fd);

    eval_unpolarized_feature_value(feature, rho, sigma + h_sigma, tau, &plus);
    eval_unpolarized_feature_value(feature, rho, sigma - h_sigma, tau, &minus);
    fd = (plus - minus) / (2.0 * h_sigma);
    if(fabs(dsigma[feature] - fd) > 1e-8 * (fabs(fd) > 1.0 ? fabs(fd) : 1.0)){
      fprintf(stderr, "finite difference sigma feature %zu analytic %.17g fd %.17g plus %.17g minus %.17g\n",
              feature, dsigma[feature], fd, plus, minus);
    }
    check_close(dsigma[feature], fd);

    eval_unpolarized_feature_value(feature, rho, sigma, tau + h_tau, &plus);
    eval_unpolarized_feature_value(feature, rho, sigma, tau - h_tau, &minus);
    fd = (plus - minus) / (2.0 * h_tau);
    if(fabs(dtau[feature] - fd) > 1e-8 * (fabs(fd) > 1.0 ? fabs(fd) : 1.0)){
      fprintf(stderr, "finite difference tau feature %zu analytic %.17g fd %.17g plus %.17g minus %.17g\n",
              feature, dtau[feature], fd, plus, minus);
    }
    check_close(dtau[feature], fd);
  }
}

int
main(void)
{
  check_metadata();
  check_dimensions();
  check_handle();
  check_error_paths();
  check_output_size_matrix();
  check_unpolarized_eval();
  check_unpolarized_order0_only_eval();
  check_polarized_eval();
  check_polarized_ml25_sample_values();
  check_deterministic_grid_for_spin(IFXC_UNPOLARIZED);
  check_deterministic_grid_for_spin(IFXC_POLARIZED);
  check_finite_difference_first_derivatives();
  return 0;
}
