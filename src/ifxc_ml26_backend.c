#include "ifxc_internal.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "ifxc_mgga.h"

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
           IFXC_MGGA_FLAGS_HAVE_EXC,
  .dens_threshold = 1e-15
};

#define IFXC_ML26_AD_DIM 7

typedef struct {
  double value;
  double gradient[IFXC_ML26_AD_DIM];
  double hessian[IFXC_ML26_AD_DIM][IFXC_ML26_AD_DIM];
} ifxc_ml26_ad2;

static ifxc_ml26_ad2
ifxc_ml26_ad_constant(double value)
{
  ifxc_ml26_ad2 result = {0};
  result.value = value;
  return result;
}

static ifxc_ml26_ad2
ifxc_ml26_ad_variable(double value, size_t variable)
{
  ifxc_ml26_ad2 result = {0};
  result.value = value;
  result.gradient[variable] = 1.0;
  return result;
}

static ifxc_ml26_ad2
ifxc_ml26_ad_add(ifxc_ml26_ad2 left, ifxc_ml26_ad2 right)
{
  ifxc_ml26_ad2 result;
  size_t i;
  size_t j;

  result.value = left.value + right.value;
  for(i = 0; i < IFXC_ML26_AD_DIM; ++i){
    result.gradient[i] = left.gradient[i] + right.gradient[i];
    for(j = 0; j < IFXC_ML26_AD_DIM; ++j){
      result.hessian[i][j] =
          left.hessian[i][j] + right.hessian[i][j];
    }
  }
  return result;
}

static ifxc_ml26_ad2
ifxc_ml26_ad_scale(ifxc_ml26_ad2 value, double scale)
{
  ifxc_ml26_ad2 result;
  size_t i;
  size_t j;

  result.value = scale * value.value;
  for(i = 0; i < IFXC_ML26_AD_DIM; ++i){
    result.gradient[i] = scale * value.gradient[i];
    for(j = 0; j < IFXC_ML26_AD_DIM; ++j){
      result.hessian[i][j] = scale * value.hessian[i][j];
    }
  }
  return result;
}

static ifxc_ml26_ad2
ifxc_ml26_ad_multiply(ifxc_ml26_ad2 left, ifxc_ml26_ad2 right)
{
  ifxc_ml26_ad2 result;
  size_t i;
  size_t j;

  result.value = left.value * right.value;
  for(i = 0; i < IFXC_ML26_AD_DIM; ++i){
    result.gradient[i] =
        left.gradient[i] * right.value +
        left.value * right.gradient[i];
    for(j = 0; j < IFXC_ML26_AD_DIM; ++j){
      result.hessian[i][j] =
          left.hessian[i][j] * right.value +
          left.gradient[i] * right.gradient[j] +
          left.gradient[j] * right.gradient[i] +
          left.value * right.hessian[i][j];
    }
  }
  return result;
}

static ifxc_ml26_ad2
ifxc_ml26_ad_unary(
    ifxc_ml26_ad2 argument,
    double value,
    double first,
    double second)
{
  ifxc_ml26_ad2 result;
  size_t i;
  size_t j;

  result.value = value;
  for(i = 0; i < IFXC_ML26_AD_DIM; ++i){
    result.gradient[i] = first * argument.gradient[i];
    for(j = 0; j < IFXC_ML26_AD_DIM; ++j){
      result.hessian[i][j] =
          second * argument.gradient[i] * argument.gradient[j] +
          first * argument.hessian[i][j];
    }
  }
  return result;
}

static ifxc_ml26_ad2
ifxc_ml26_ad_power(ifxc_ml26_ad2 argument, double exponent)
{
  double value = pow(argument.value, exponent);
  double first = exponent * pow(argument.value, exponent - 1.0);
  double second =
      exponent * (exponent - 1.0) *
      pow(argument.value, exponent - 2.0);
  return ifxc_ml26_ad_unary(argument, value, first, second);
}

static ifxc_ml26_ad2
ifxc_ml26_ad_inverse(ifxc_ml26_ad2 argument)
{
  double inverse = 1.0 / argument.value;
  return ifxc_ml26_ad_unary(
      argument, inverse, -inverse * inverse,
      2.0 * inverse * inverse * inverse);
}

static ifxc_ml26_ad2
ifxc_ml26_ad_divide(ifxc_ml26_ad2 numerator, ifxc_ml26_ad2 denominator)
{
  return ifxc_ml26_ad_multiply(
      numerator, ifxc_ml26_ad_inverse(denominator));
}

static ifxc_ml26_ad2
ifxc_ml26_ad_floor(ifxc_ml26_ad2 value, double floor_value)
{
  if(value.value < floor_value){
    return ifxc_ml26_ad_constant(floor_value);
  }
  return value;
}

static ifxc_ml26_ad2
ifxc_ml26_cs1_density_factor(ifxc_ml26_ad2 density)
{
  const double d = 0.349;
  ifxc_ml26_ad2 root = ifxc_ml26_ad_power(density, 1.0 / 3.0);
  return ifxc_ml26_ad_divide(
      root, ifxc_ml26_ad_add(root, ifxc_ml26_ad_constant(d)));
}

static ifxc_ml26_ad2
ifxc_ml26_cs1_gradient_factor(
    ifxc_ml26_ad2 density,
    ifxc_ml26_ad2 sigma)
{
  const double gamma = 0.006;
  ifxc_ml26_ad2 sigma_squared = ifxc_ml26_ad_multiply(sigma, sigma);
  ifxc_ml26_ad2 denominator = ifxc_ml26_ad_add(
      ifxc_ml26_ad_power(density, 8.0 / 3.0),
      ifxc_ml26_ad_scale(sigma, gamma));
  return ifxc_ml26_ad_scale(
      ifxc_ml26_ad_divide(
          sigma_squared,
          ifxc_ml26_ad_multiply(denominator, denominator)),
      gamma * gamma);
}

static ifxc_ml26_ad2
ifxc_ml26_cs1_spin_fraction(
    ifxc_ml26_ad2 spin_density,
    ifxc_ml26_ad2 total_density)
{
  ifxc_ml26_ad2 fraction =
      ifxc_ml26_ad_divide(spin_density, total_density);
  if(2.0 * fraction.value <= DBL_EPSILON){
    return ifxc_ml26_ad_constant(0.5 * DBL_EPSILON);
  }
  return fraction;
}

static void
ifxc_ml26_eval_cs1_point(
    const ifxc_handle_impl *impl,
    const ifxc_input *input,
    size_t point,
    ifxc_ml26_ad2 features[3])
{
  const double c1 = -0.018897;
  const double c2 = 0.155240;
  const double c3 = -0.159068;
  const double c4 = 0.007953;
  const double sigma_floor = pow(1e-15, 8.0 / 3.0);
  ifxc_ml26_ad2 raw_total;
  ifxc_ml26_ad2 density_a;
  ifxc_ml26_ad2 density_b;
  ifxc_ml26_ad2 sigma_aa;
  ifxc_ml26_ad2 sigma_ab;
  ifxc_ml26_ad2 sigma_bb;
  ifxc_ml26_ad2 total_density;
  ifxc_ml26_ad2 total_sigma;
  ifxc_ml26_ad2 fraction_a;
  ifxc_ml26_ad2 fraction_b;
  ifxc_ml26_ad2 density_factor_a;
  ifxc_ml26_ad2 density_factor_b;
  ifxc_ml26_ad2 same_spin_1;
  ifxc_ml26_ad2 same_spin_2;
  ifxc_ml26_ad2 opposite_spin;
  ifxc_ml26_ad2 opposite_prefactor;
  ifxc_ml26_ad2 opposite_bracket;

  features[0] = ifxc_ml26_ad_constant(0.0);
  features[1] = ifxc_ml26_ad_constant(0.0);
  features[2] = ifxc_ml26_ad_constant(0.0);

  if(impl->nspin == IFXC_UNPOLARIZED){
    ifxc_ml26_ad2 density =
        ifxc_ml26_ad_variable(input->rho[point], 0);
    ifxc_ml26_ad2 sigma =
        ifxc_ml26_ad_variable(input->sigma[point], 1);
    if(density.value < 1e-15){
      return;
    }
    raw_total = density;
    density = ifxc_ml26_ad_floor(density, 1e-15);
    sigma = ifxc_ml26_ad_floor(sigma, sigma_floor);
    density_a = ifxc_ml26_ad_scale(density, 0.5);
    density_b = density_a;
    sigma_aa = ifxc_ml26_ad_scale(sigma, 0.25);
    sigma_ab = sigma_aa;
    sigma_bb = sigma_aa;
  }else{
    ifxc_ml26_ad2 raw_density_a =
        ifxc_ml26_ad_variable(input->rho[point], 0);
    ifxc_ml26_ad2 raw_density_b =
        ifxc_ml26_ad_variable(
            input->rho[input->npoints + point], 1);
    ifxc_ml26_ad2 raw_sigma_aa =
        ifxc_ml26_ad_variable(input->sigma[point], 2);
    ifxc_ml26_ad2 raw_sigma_ab =
        ifxc_ml26_ad_variable(
            input->sigma[input->npoints + point], 3);
    ifxc_ml26_ad2 raw_sigma_bb =
        ifxc_ml26_ad_variable(
            input->sigma[2 * input->npoints + point], 4);
    ifxc_ml26_ad2 sigma_average;

    raw_total = ifxc_ml26_ad_add(raw_density_a, raw_density_b);
    if(raw_total.value < 1e-15){
      return;
    }
    density_a = ifxc_ml26_ad_floor(raw_density_a, 1e-15);
    density_b = ifxc_ml26_ad_floor(raw_density_b, 1e-15);
    sigma_aa = ifxc_ml26_ad_floor(raw_sigma_aa, sigma_floor);
    sigma_bb = ifxc_ml26_ad_floor(raw_sigma_bb, sigma_floor);
    sigma_average = ifxc_ml26_ad_scale(
        ifxc_ml26_ad_add(sigma_aa, sigma_bb), 0.5);
    if(raw_sigma_ab.value < -sigma_average.value){
      sigma_ab = ifxc_ml26_ad_scale(sigma_average, -1.0);
    }else if(raw_sigma_ab.value > sigma_average.value){
      sigma_ab = sigma_average;
    }else{
      sigma_ab = raw_sigma_ab;
    }
  }

  total_density = ifxc_ml26_ad_add(density_a, density_b);
  total_sigma = ifxc_ml26_ad_add(
      ifxc_ml26_ad_add(sigma_aa, sigma_bb),
      ifxc_ml26_ad_scale(sigma_ab, 2.0));
  fraction_a =
      ifxc_ml26_cs1_spin_fraction(density_a, total_density);
  fraction_b =
      ifxc_ml26_cs1_spin_fraction(density_b, total_density);
  density_factor_a = ifxc_ml26_cs1_density_factor(density_a);
  density_factor_b = ifxc_ml26_cs1_density_factor(density_b);

  same_spin_1 = ifxc_ml26_ad_add(
      ifxc_ml26_ad_multiply(fraction_a, density_factor_a),
      ifxc_ml26_ad_multiply(fraction_b, density_factor_b));
  features[0] = ifxc_ml26_ad_scale(
      ifxc_ml26_ad_multiply(raw_total, same_spin_1), c1);

  same_spin_2 = ifxc_ml26_ad_add(
      ifxc_ml26_ad_multiply(
          ifxc_ml26_ad_multiply(fraction_a, density_factor_a),
          ifxc_ml26_cs1_gradient_factor(density_a, sigma_aa)),
      ifxc_ml26_ad_multiply(
          ifxc_ml26_ad_multiply(fraction_b, density_factor_b),
          ifxc_ml26_cs1_gradient_factor(density_b, sigma_bb)));
  features[1] = ifxc_ml26_ad_scale(
      ifxc_ml26_ad_multiply(raw_total, same_spin_2), c2);

  opposite_prefactor = ifxc_ml26_ad_multiply(
      ifxc_ml26_ad_divide(
          ifxc_ml26_ad_multiply(density_a, density_b),
          ifxc_ml26_ad_multiply(total_density, total_density)),
      ifxc_ml26_cs1_density_factor(total_density));
  opposite_bracket = ifxc_ml26_ad_add(
      ifxc_ml26_ad_constant(c3),
      ifxc_ml26_ad_scale(
          ifxc_ml26_cs1_gradient_factor(total_density, total_sigma),
          c4));
  opposite_spin = ifxc_ml26_ad_multiply(
      opposite_prefactor, opposite_bracket);
  features[2] = ifxc_ml26_ad_multiply(raw_total, opposite_spin);
}

static size_t
ifxc_ml26_variable_dimension(
    const ifxc_dimensions_t *dims,
    ifxc_variable variable)
{
  switch(variable){
  case IFXC_VAR_RHO:
    return dims->rho;
  case IFXC_VAR_SIGMA:
    return dims->sigma;
  case IFXC_VAR_TAU:
    return dims->tau;
  default:
    return 0;
  }
}

static size_t
ifxc_ml26_variable_offset(
    const ifxc_dimensions_t *dims,
    ifxc_variable variable)
{
  switch(variable){
  case IFXC_VAR_RHO:
    return 0;
  case IFXC_VAR_SIGMA:
    return dims->rho;
  case IFXC_VAR_TAU:
    return dims->rho + dims->sigma;
  default:
    return IFXC_ML26_AD_DIM;
  }
}

static double
ifxc_ml26_cs1_derivative_component(
    const ifxc_dimensions_t *dims,
    const ifxc_deriv_entry *entry,
    const ifxc_ml26_ad2 *feature,
    size_t component)
{
  size_t first_dimension;
  size_t first_offset;

  first_dimension =
      ifxc_ml26_variable_dimension(dims, entry->vars[0]);
  first_offset = ifxc_ml26_variable_offset(dims, entry->vars[0]);
  if(entry->order == 1){
    return feature->gradient[first_offset + component];
  }

  if(entry->vars[0] == entry->vars[1]){
    size_t first;
    size_t second;
    size_t current = 0;
    for(first = 0; first < first_dimension; ++first){
      for(second = first; second < first_dimension; ++second){
        if(current == component){
          return feature->hessian[
              first_offset + first][first_offset + second];
        }
        ++current;
      }
    }
    return 0.0;
  }

  {
    size_t second_dimension =
        ifxc_ml26_variable_dimension(dims, entry->vars[1]);
    size_t second_offset =
        ifxc_ml26_variable_offset(dims, entry->vars[1]);
    size_t first = component / second_dimension;
    size_t second = component % second_dimension;
    return feature->hessian[
        first_offset + first][second_offset + second];
  }
}

static ifxc_status
ifxc_ml26_eval_derivative_entry(
    const ifxc_handle_impl *impl,
    const ifxc_input *input,
    const ifxc_deriv_entry *entry)
{
  ifxc_handle_impl prefix_impl = *impl;
  ifxc_func_type prefix_func;
  ifxc_deriv_entry prefix_entry = *entry;
  double *prefix_output = NULL;
  size_t ncomponents;
  size_t point;
  size_t feature;
  ifxc_status status;

  status = ifxc_derivative_component_count(
      &impl->dims, entry, &ncomponents);
  if(status != IFXC_OK){
    return status;
  }
  if(input->npoints == 0){
    return IFXC_OK;
  }
  if(input->npoints > SIZE_MAX /
      (IFXC_ML25_NFEATURES * ncomponents * sizeof(double))){
    return IFXC_E_ALLOCATION;
  }
  prefix_output = (double *)malloc(
      input->npoints * IFXC_ML25_NFEATURES *
      ncomponents * sizeof(double));
  if(prefix_output == NULL){
    return IFXC_E_ALLOCATION;
  }

  prefix_impl.feature_set = IFXC_FEATURE_SET_ML25;
  prefix_impl.nfeatures = IFXC_ML25_NFEATURES;
  prefix_impl.max_deriv_order = IFXC_ML25_GENERATED_MAX_ORDER;
  prefix_func.impl = &prefix_impl;
  prefix_entry.out = prefix_output;
  status = ifxc_ml25_eval(
      &prefix_func, input, 1, &prefix_entry);
  if(status != IFXC_OK){
    free(prefix_output);
    return status;
  }

  for(point = 0; point < input->npoints; ++point){
    ifxc_ml26_ad2 cs1[3];
    size_t component;
    double weight =
        (entry->target == IFXC_TARGET_INTEGRAL)
        ? input->weights[point]
        : 1.0;
    ifxc_ml26_eval_cs1_point(impl, input, point, cs1);
    for(feature = 0; feature < IFXC_ML25_NFEATURES; ++feature){
      for(component = 0; component < ncomponents; ++component){
        entry->out[
            (point * IFXC_ML26_NFEATURES + feature) *
                ncomponents + component] =
            prefix_output[
                (point * IFXC_ML25_NFEATURES + feature) *
                    ncomponents + component];
      }
    }
    for(feature = 0; feature < 3; ++feature){
      for(component = 0; component < ncomponents; ++component){
        entry->out[
            (point * IFXC_ML26_NFEATURES +
             IFXC_ML25_NFEATURES + feature) *
                ncomponents + component] =
            weight * ifxc_ml26_cs1_derivative_component(
                &impl->dims, entry, &cs1[feature], component);
      }
    }
  }

  free(prefix_output);
  return IFXC_OK;
}

static void
ifxc_ml26_prepare_private_func(
    ifxc_mgga_func_type *p,
    ifxc_nspin nspin)
{
  memset(p, 0, sizeof(*p));
  p->info = &IFXC_ML26_PRIVATE_INFO;
  p->nspin = (nspin == IFXC_POLARIZED)
      ? IFXC_MGGA_POLARIZED
      : IFXC_MGGA_UNPOLARIZED;
  p->dim.zk = IFXC_ML26_NFEATURES;
  if(p->nspin == IFXC_MGGA_POLARIZED){
    p->dim.rho = 2;
    p->dim.sigma = 3;
    p->dim.tau = 2;
  }else{
    p->dim.rho = 1;
    p->dim.sigma = 1;
    p->dim.tau = 1;
  }
  p->dens_threshold = p->info->dens_threshold;
  p->sigma_threshold = pow(p->dens_threshold, 4.0 / 3.0);
  p->zeta_threshold = DBL_EPSILON;
  p->tau_threshold = 1e-20;
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
  ifxc_mgga_func_type p;
  size_t entry_index;

  if(impl == NULL || input == NULL || entries == NULL){
    return IFXC_E_INVALID_ARGUMENT;
  }
  if(impl->feature_set != IFXC_FEATURE_SET_ML26){
    return IFXC_E_UNKNOWN_FEATURE_SET;
  }

  ifxc_ml26_prepare_private_func(&p, impl->nspin);
  for(entry_index = 0; entry_index < nentries; ++entry_index){
    const ifxc_deriv_entry *entry = &entries[entry_index];
    if(entry->order > 0){
      ifxc_status status =
          ifxc_ml26_eval_derivative_entry(impl, input, entry);
      if(status != IFXC_OK){
        return status;
      }
    }else if(entry->target == IFXC_TARGET_LOCAL){
      ifxc_status status = ifxc_ml26_combined_eval_order0_component_major(
          &p, input->npoints, input->rho, input->sigma, input->tau,
          NULL, entry->out, NULL);
      if(status != IFXC_OK){
        return status;
      }
    }else if(entry->target == IFXC_TARGET_INTEGRAL){
      ifxc_status status = ifxc_ml26_combined_eval_order0_component_major(
          &p, input->npoints, input->rho, input->sigma, input->tau,
          input->weights, NULL, entry->out);
      if(status != IFXC_OK){
        return status;
      }
    }else{
      return IFXC_E_INVALID_ARGUMENT;
    }
  }

  return IFXC_OK;
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
  ifxc_handle_impl prefix_impl;
  size_t point;
  size_t component;
  ifxc_status status;

  if(impl == NULL || input == NULL || coeffs == NULL ||
     d_rho == NULL || d_sigma == NULL || d_tau == NULL){
    return IFXC_E_INVALID_ARGUMENT;
  }
  if(impl->feature_set != IFXC_FEATURE_SET_ML26){
    return IFXC_E_UNKNOWN_FEATURE_SET;
  }

  prefix_impl = *impl;
  prefix_impl.feature_set = IFXC_FEATURE_SET_ML25;
  prefix_impl.nfeatures = IFXC_ML25_NFEATURES;
  prefix_impl.max_deriv_order = IFXC_ML25_GENERATED_MAX_ORDER;
  status = ifxc_ml25_eval_first_derivatives_contracted(
      &prefix_impl, input, coeffs, d_rho, d_sigma, d_tau);
  if(status != IFXC_OK){
    return status;
  }

  for(point = 0; point < input->npoints; ++point){
    ifxc_ml26_ad2 cs1[3];
    ifxc_ml26_eval_cs1_point(impl, input, point, cs1);
    for(component = 0; component < impl->dims.rho; ++component){
      size_t variable = component;
      d_rho[component * input->npoints + point] +=
          coeffs[IFXC_ML26_CS1_SAME_SPIN_1] *
              cs1[0].gradient[variable] +
          coeffs[IFXC_ML26_CS1_SAME_SPIN_2] *
              cs1[1].gradient[variable] +
          coeffs[IFXC_ML26_CS1_OPPOSITE_SPIN] *
              cs1[2].gradient[variable];
    }
    for(component = 0; component < impl->dims.sigma; ++component){
      size_t variable = impl->dims.rho + component;
      d_sigma[component * input->npoints + point] +=
          coeffs[IFXC_ML26_CS1_SAME_SPIN_1] *
              cs1[0].gradient[variable] +
          coeffs[IFXC_ML26_CS1_SAME_SPIN_2] *
              cs1[1].gradient[variable] +
          coeffs[IFXC_ML26_CS1_OPPOSITE_SPIN] *
              cs1[2].gradient[variable];
    }
  }

  return IFXC_OK;
}
