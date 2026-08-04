/* SPDX-License-Identifier: MPL-2.0 */

#include "ifxc_internal.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "ifxc_mgga.h"
#include "ifxc_vector_mgga.h"

ifxc_mgga_funcs_variants *ifxc_ml25_combined_work_mgga(void);
int ifxc_ml25_has_combined_order(unsigned int order);
ifxc_status ifxc_ml25_combined_eval_contracted_component_major(
    const ifxc_mgga_func_type *p,
    size_t npoints,
    const double *rho,
    const double *sigma,
    const double *tau,
    double *d_rho,
    double *d_sigma,
    double *d_tau);
ifxc_status ifxc_ml25_combined_integral_order0_component_major(
    const ifxc_mgga_func_type *p,
    size_t npoints,
    const double *rho,
    const double *sigma,
    const double *tau,
    const double *weights,
    double *features);

static const ifxc_mgga_func_info IFXC_ML25_PRIVATE_INFO = {
  .flags = IFXC_MGGA_FLAGS_3D | IFXC_MGGA_FLAGS_NEEDS_TAU | IFXC_MGGA_FLAGS_HAVE_ALL,
  .dens_threshold = 1e-15
};

#define IFXC_ML25_PRIVATE_FIELDS(X) \
  X(zk) \
  X(vrho) \
  X(vsigma) \
  X(vtau) \
  X(v2rho2) \
  X(v2rhosigma) \
  X(v2rhotau) \
  X(v2sigma2) \
  X(v2sigmatau) \
  X(v2tau2) \
  X(v3rho3) \
  X(v3rho2sigma) \
  X(v3rho2tau) \
  X(v3rhosigma2) \
  X(v3rhosigmatau) \
  X(v3rhotau2) \
  X(v3sigma3) \
  X(v3sigma2tau) \
  X(v3sigmatau2) \
  X(v3tau3) \
  X(v4rho4) \
  X(v4rho3sigma) \
  X(v4rho3tau) \
  X(v4rho2sigma2) \
  X(v4rho2sigmatau) \
  X(v4rho2tau2) \
  X(v4rhosigma3) \
  X(v4rhosigma2tau) \
  X(v4rhosigmatau2) \
  X(v4rhotau3) \
  X(v4sigma4) \
  X(v4sigma3tau) \
  X(v4sigma2tau2) \
  X(v4sigmatau3) \
  X(v4tau4)

void ifxc_vector_mgga_free_inputs(double *rho_tm, double *sigma_tm, double *lapl_tm, double *tau_tm);
void ifxc_vector_mgga_free_private_outputs(ifxc_mgga_out_params *out);
static ifxc_status ifxc_ml25_eval_combined(
    const ifxc_handle_impl *impl,
    const ifxc_input *input,
    size_t nentries,
    const ifxc_deriv_entry *entries,
    unsigned int max_order);
static double
ifxc_ml25_total_density_at_point(const ifxc_input *input,
                                 const ifxc_dimensions_t *dims,
                                 size_t point)
{
  double dens = 0.0;
  size_t comp;

  for(comp = 0; comp < dims->rho; ++comp){
    dens += input->rho[comp * input->npoints + point];
  }

  return dens;
}

static double
ifxc_ml25_order0_kernel(
    const ifxc_input *input,
    const ifxc_dimensions_t *dims,
    size_t point,
    double feature_scale,
    double generated_zk_value)
{
  return feature_scale *
         ifxc_ml25_total_density_at_point(input, dims, point) *
         generated_zk_value;
}

void
ifxc_vector_mgga_set_private_dimensions(int nspin, ifxc_mgga_dimensions *dim)
{
  memset(dim, 0, sizeof(*dim));
  dim->zk = 1;

  if(nspin == IFXC_MGGA_UNPOLARIZED){
    dim->rho = dim->sigma = dim->tau = 1;
    dim->vrho = dim->vsigma = dim->vtau = 1;
    dim->v2rho2 = dim->v2rhosigma = dim->v2rhotau = 1;
    dim->v2sigma2 = dim->v2sigmatau = dim->v2tau2 = 1;
    dim->v3rho3 = dim->v3rho2sigma = dim->v3rho2tau = 1;
    dim->v3rhosigma2 = dim->v3rhosigmatau = dim->v3rhotau2 = 1;
    dim->v3sigma3 = dim->v3sigma2tau = dim->v3sigmatau2 = 1;
    dim->v3tau3 = 1;
    dim->v4rho4 = dim->v4rho3sigma = dim->v4rho3tau = 1;
    dim->v4rho2sigma2 = dim->v4rho2sigmatau = dim->v4rho2tau2 = 1;
    dim->v4rhosigma3 = dim->v4rhosigma2tau = dim->v4rhosigmatau2 = 1;
    dim->v4rhotau3 = dim->v4sigma4 = dim->v4sigma3tau = 1;
    dim->v4sigma2tau2 = dim->v4sigmatau3 = dim->v4tau4 = 1;
    return;
  }

  dim->rho = 2;
  dim->sigma = 3;
  dim->tau = 2;
  dim->vrho = 2;
  dim->vsigma = 3;
  dim->vtau = 2;

  dim->v2rho2 = 3;
  dim->v2rhosigma = 2 * 3;
  dim->v2rhotau = 2 * 2;
  dim->v2sigma2 = 6;
  dim->v2sigmatau = 3 * 2;
  dim->v2tau2 = 3;

  dim->v3rho3 = 4;
  dim->v3rho2sigma = 3 * 3;
  dim->v3rho2tau = 3 * 2;
  dim->v3rhosigma2 = 2 * 6;
  dim->v3rhosigmatau = 2 * 3 * 2;
  dim->v3rhotau2 = 2 * 3;
  dim->v3sigma3 = 10;
  dim->v3sigma2tau = 6 * 2;
  dim->v3sigmatau2 = 3 * 3;
  dim->v3tau3 = 4;

  dim->v4rho4 = 5;
  dim->v4rho3sigma = 4 * 3;
  dim->v4rho3tau = 4 * 2;
  dim->v4rho2sigma2 = 3 * 6;
  dim->v4rho2sigmatau = 3 * 3 * 2;
  dim->v4rho2tau2 = 3 * 3;
  dim->v4rhosigma3 = 2 * 10;
  dim->v4rhosigma2tau = 2 * 6 * 2;
  dim->v4rhosigmatau2 = 2 * 3 * 3;
  dim->v4rhotau3 = 2 * 4;
  dim->v4sigma4 = 15;
  dim->v4sigma3tau = 10 * 2;
  dim->v4sigma2tau2 = 6 * 3;
  dim->v4sigmatau3 = 3 * 4;
  dim->v4tau4 = 5;
}

const double *
ifxc_vector_mgga_select_source_buffer(const ifxc_mgga_out_params *out, const ifxc_deriv_entry *entry)
{
  switch(entry->order){
  case 0:
    return out->zk;

  case 1:
    switch(entry->vars[0]){
    case IFXC_VAR_RHO: return out->vrho;
    case IFXC_VAR_SIGMA: return out->vsigma;
    case IFXC_VAR_TAU: return out->vtau;
    default: return NULL;
    }

  case 2:
    switch(entry->vars[0]){
    case IFXC_VAR_RHO:
      switch(entry->vars[1]){
      case IFXC_VAR_RHO: return out->v2rho2;
      case IFXC_VAR_SIGMA: return out->v2rhosigma;
      case IFXC_VAR_TAU: return out->v2rhotau;
      default: return NULL;
      }
    case IFXC_VAR_SIGMA:
      switch(entry->vars[1]){
      case IFXC_VAR_SIGMA: return out->v2sigma2;
      case IFXC_VAR_TAU: return out->v2sigmatau;
      default: return NULL;
      }
    case IFXC_VAR_TAU:
      return (entry->vars[1] == IFXC_VAR_TAU) ? out->v2tau2 : NULL;
    default:
      return NULL;
    }

  case 3:
    switch(entry->vars[0]){
    case IFXC_VAR_RHO:
      switch(entry->vars[1]){
      case IFXC_VAR_RHO:
        switch(entry->vars[2]){
        case IFXC_VAR_RHO: return out->v3rho3;
        case IFXC_VAR_SIGMA: return out->v3rho2sigma;
        case IFXC_VAR_TAU: return out->v3rho2tau;
        default: return NULL;
        }
      case IFXC_VAR_SIGMA:
        switch(entry->vars[2]){
        case IFXC_VAR_SIGMA: return out->v3rhosigma2;
        case IFXC_VAR_TAU: return out->v3rhosigmatau;
        default: return NULL;
        }
      case IFXC_VAR_TAU:
        return (entry->vars[2] == IFXC_VAR_TAU) ? out->v3rhotau2 : NULL;
      default:
        return NULL;
      }
    case IFXC_VAR_SIGMA:
      switch(entry->vars[1]){
      case IFXC_VAR_SIGMA:
        switch(entry->vars[2]){
        case IFXC_VAR_SIGMA: return out->v3sigma3;
        case IFXC_VAR_TAU: return out->v3sigma2tau;
        default: return NULL;
        }
      case IFXC_VAR_TAU:
        return (entry->vars[2] == IFXC_VAR_TAU) ? out->v3sigmatau2 : NULL;
      default:
        return NULL;
      }
    case IFXC_VAR_TAU:
      return (entry->vars[1] == IFXC_VAR_TAU && entry->vars[2] == IFXC_VAR_TAU)
        ? out->v3tau3
        : NULL;
    default:
      return NULL;
    }

  case 4:
    switch(entry->vars[0]){
    case IFXC_VAR_RHO:
      switch(entry->vars[1]){
      case IFXC_VAR_RHO:
        switch(entry->vars[2]){
        case IFXC_VAR_RHO:
          switch(entry->vars[3]){
          case IFXC_VAR_RHO: return out->v4rho4;
          case IFXC_VAR_SIGMA: return out->v4rho3sigma;
          case IFXC_VAR_TAU: return out->v4rho3tau;
          default: return NULL;
          }
        case IFXC_VAR_SIGMA:
          switch(entry->vars[3]){
          case IFXC_VAR_SIGMA: return out->v4rho2sigma2;
          case IFXC_VAR_TAU: return out->v4rho2sigmatau;
          default: return NULL;
          }
        case IFXC_VAR_TAU:
          return (entry->vars[3] == IFXC_VAR_TAU) ? out->v4rho2tau2 : NULL;
        default:
          return NULL;
        }
      case IFXC_VAR_SIGMA:
        switch(entry->vars[2]){
        case IFXC_VAR_SIGMA:
          switch(entry->vars[3]){
          case IFXC_VAR_SIGMA: return out->v4rhosigma3;
          case IFXC_VAR_TAU: return out->v4rhosigma2tau;
          default: return NULL;
          }
        case IFXC_VAR_TAU:
          return (entry->vars[3] == IFXC_VAR_TAU) ? out->v4rhosigmatau2 : NULL;
        default:
          return NULL;
        }
      case IFXC_VAR_TAU:
        return (entry->vars[2] == IFXC_VAR_TAU && entry->vars[3] == IFXC_VAR_TAU)
          ? out->v4rhotau3
          : NULL;
      default:
        return NULL;
      }

    case IFXC_VAR_SIGMA:
      switch(entry->vars[1]){
      case IFXC_VAR_SIGMA:
        switch(entry->vars[2]){
        case IFXC_VAR_SIGMA:
          switch(entry->vars[3]){
          case IFXC_VAR_SIGMA: return out->v4sigma4;
          case IFXC_VAR_TAU: return out->v4sigma3tau;
          default: return NULL;
          }
        case IFXC_VAR_TAU:
          return (entry->vars[3] == IFXC_VAR_TAU) ? out->v4sigma2tau2 : NULL;
        default:
          return NULL;
        }
      case IFXC_VAR_TAU:
        return (entry->vars[2] == IFXC_VAR_TAU && entry->vars[3] == IFXC_VAR_TAU)
          ? out->v4sigmatau3
          : NULL;
      default:
        return NULL;
      }

    case IFXC_VAR_TAU:
      return (entry->vars[1] == IFXC_VAR_TAU &&
              entry->vars[2] == IFXC_VAR_TAU &&
              entry->vars[3] == IFXC_VAR_TAU)
        ? out->v4tau4
        : NULL;

    default:
      return NULL;
    }

  default:
    return NULL;
  }
}

void
ifxc_vector_mgga_zero_dimensions_above_order(ifxc_mgga_dimensions *dims, unsigned int max_order)
{
  if(max_order < 4){
    dims->v4rho4 = dims->v4rho3sigma = dims->v4rho3lapl = dims->v4rho3tau = 0;
    dims->v4rho2sigma2 = dims->v4rho2sigmalapl = dims->v4rho2sigmatau = 0;
    dims->v4rho2lapl2 = dims->v4rho2lapltau = dims->v4rho2tau2 = 0;
    dims->v4rhosigma3 = dims->v4rhosigma2lapl = dims->v4rhosigma2tau = 0;
    dims->v4rhosigmalapl2 = dims->v4rhosigmalapltau = dims->v4rhosigmatau2 = 0;
    dims->v4rholapl3 = dims->v4rholapl2tau = dims->v4rholapltau2 = 0;
    dims->v4rhotau3 = dims->v4sigma4 = dims->v4sigma3lapl = 0;
    dims->v4sigma3tau = dims->v4sigma2lapl2 = dims->v4sigma2lapltau = 0;
    dims->v4sigma2tau2 = dims->v4sigmalapl3 = dims->v4sigmalapl2tau = 0;
    dims->v4sigmalapltau2 = dims->v4sigmatau3 = dims->v4lapl4 = 0;
    dims->v4lapl3tau = dims->v4lapl2tau2 = dims->v4lapltau3 = dims->v4tau4 = 0;
  }
  if(max_order < 3){
    dims->v3rho3 = dims->v3rho2sigma = dims->v3rho2lapl = dims->v3rho2tau = 0;
    dims->v3rhosigma2 = dims->v3rhosigmalapl = dims->v3rhosigmatau = 0;
    dims->v3rholapl2 = dims->v3rholapltau = dims->v3rhotau2 = 0;
    dims->v3sigma3 = dims->v3sigma2lapl = dims->v3sigma2tau = 0;
    dims->v3sigmalapl2 = dims->v3sigmalapltau = dims->v3sigmatau2 = 0;
    dims->v3lapl3 = dims->v3lapl2tau = dims->v3lapltau2 = dims->v3tau3 = 0;
  }
  if(max_order < 2){
    dims->v2rho2 = dims->v2rhosigma = dims->v2rholapl = dims->v2rhotau = 0;
    dims->v2sigma2 = dims->v2sigmalapl = dims->v2sigmatau = 0;
    dims->v2lapl2 = dims->v2lapltau = dims->v2tau2 = 0;
  }
  if(max_order < 1){
    dims->vrho = dims->vsigma = dims->vlapl = dims->vtau = 0;
  }
}

void
ifxc_vector_mgga_scale_output_dimensions(ifxc_mgga_dimensions *dims, size_t nfeatures)
{
#define IFXC_SCALE_FIELD(name) dims->name *= (int)nfeatures;

  IFXC_ML25_PRIVATE_FIELDS(IFXC_SCALE_FIELD)

#undef IFXC_SCALE_FIELD
}

ifxc_status
ifxc_vector_mgga_copy_combined_entry(
    const ifxc_input *input,
    const ifxc_dimensions_t *dims,
    size_t nfeatures,
    const ifxc_deriv_entry *entry,
    const double *src,
    double scale)
{
  size_t npoints = input->npoints;
  size_t ncomp = 0;
  size_t point;
  size_t feature;
  ifxc_status status = ifxc_derivative_component_count(dims, entry, &ncomp);
  if(status != IFXC_OK){
    return status;
  }

  if(entry->target == IFXC_TARGET_LOCAL){
    for(point = 0; point < npoints; ++point){
      for(feature = 0; feature < nfeatures; ++feature){
        size_t comp;
        for(comp = 0; comp < ncomp; ++comp){
          size_t src_index = point * nfeatures * ncomp + comp * nfeatures + feature;
          size_t dest_index = (entry->order == 0)
            ? point * nfeatures + feature
            : (point * nfeatures + feature) * ncomp + comp;
          entry->out[dest_index] = (entry->order == 0)
            ? ifxc_ml25_order0_kernel(input, dims, point, scale, src[src_index])
            : scale * src[src_index];
        }
      }
    }
    return IFXC_OK;
  }

  if(entry->target == IFXC_TARGET_INTEGRAL){
    if(entry->order == 0){
      for(feature = 0; feature < nfeatures; ++feature){
        double sum = 0.0;
        for(point = 0; point < npoints; ++point){
          sum += input->weights[point] *
                 ifxc_ml25_order0_kernel(
                     input, dims, point, scale, src[point * nfeatures + feature]);
        }
        entry->out[feature] = sum;
      }
      return IFXC_OK;
    }

    for(point = 0; point < npoints; ++point){
      for(feature = 0; feature < nfeatures; ++feature){
        size_t comp;
        for(comp = 0; comp < ncomp; ++comp){
          size_t src_index = point * nfeatures * ncomp + comp * nfeatures + feature;
          size_t dest_index = (point * nfeatures + feature) * ncomp + comp;
          entry->out[dest_index] = input->weights[point] * scale * src[src_index];
        }
      }
    }
    return IFXC_OK;
  }

  return IFXC_E_INVALID_ARGUMENT;
}

ifxc_status
ifxc_vector_mgga_transpose_inputs(
    const ifxc_dimensions_t *dims,
    const ifxc_input *input,
    double **rho_tm,
    double **sigma_tm,
    double **lapl_tm,
    double **tau_tm)
{
  size_t npoints = input->npoints;
  size_t p;
  size_t c;

  *rho_tm = NULL;
  *sigma_tm = NULL;
  *lapl_tm = NULL;
  *tau_tm = NULL;

  if(npoints == 0){
    return IFXC_OK;
  }

  *rho_tm = (double *)calloc(npoints * dims->rho, sizeof(double));
  *sigma_tm = (double *)calloc(npoints * dims->sigma, sizeof(double));
  *tau_tm = (double *)calloc(npoints * dims->tau, sizeof(double));
  *lapl_tm = (double *)calloc(1, sizeof(double));
  if(*rho_tm == NULL || *sigma_tm == NULL || *tau_tm == NULL || *lapl_tm == NULL){
    ifxc_vector_mgga_free_inputs(*rho_tm, *sigma_tm, *lapl_tm, *tau_tm);
    *rho_tm = NULL;
    *sigma_tm = NULL;
    *lapl_tm = NULL;
    *tau_tm = NULL;
    return IFXC_E_ALLOCATION;
  }

  for(c = 0; c < dims->rho; ++c){
    for(p = 0; p < npoints; ++p){
      (*rho_tm)[p * dims->rho + c] = input->rho[c * npoints + p];
    }
  }

  for(c = 0; c < dims->sigma; ++c){
    for(p = 0; p < npoints; ++p){
      (*sigma_tm)[p * dims->sigma + c] = input->sigma[c * npoints + p];
    }
  }

  for(c = 0; c < dims->tau; ++c){
    for(p = 0; p < npoints; ++p){
      (*tau_tm)[p * dims->tau + c] = input->tau[c * npoints + p];
    }
  }

  return IFXC_OK;
}

void
ifxc_vector_mgga_free_inputs(double *rho_tm, double *sigma_tm, double *lapl_tm, double *tau_tm)
{
  free(rho_tm);
  free(sigma_tm);
  free(lapl_tm);
  free(tau_tm);
}

ifxc_status
ifxc_vector_mgga_allocate_private_outputs(
    const ifxc_mgga_dimensions *dims,
    size_t npoints,
    ifxc_mgga_out_params *out)
{
  ifxc_status status = IFXC_OK;

#define IFXC_ALLOC_FIELD(name) \
  do { \
    if(npoints > 0 && dims->name > 0){ \
      out->name = (double *)calloc(npoints * (size_t)dims->name, sizeof(double)); \
      if(out->name == NULL){ \
        status = IFXC_E_ALLOCATION; \
        goto cleanup; \
      } \
    } \
  } while(0);

  IFXC_ML25_PRIVATE_FIELDS(IFXC_ALLOC_FIELD)

#undef IFXC_ALLOC_FIELD
  return IFXC_OK;

cleanup:
  ifxc_vector_mgga_free_private_outputs(out);
  memset(out, 0, sizeof(*out));
  return status;
}

void
ifxc_vector_mgga_free_private_outputs(ifxc_mgga_out_params *out)
{
#define IFXC_FREE_FIELD(name) free(out->name);

  IFXC_ML25_PRIVATE_FIELDS(IFXC_FREE_FIELD);

#undef IFXC_FREE_FIELD
}

unsigned int
ifxc_vector_mgga_max_requested_order(size_t nentries, const ifxc_deriv_entry *entries)
{
  unsigned int max_order = 0;
  size_t i;

  for(i = 0; i < nentries; ++i){
    if(entries[i].order > max_order){
      max_order = entries[i].order;
    }
  }

  return max_order;
}

static ifxc_status
ifxc_ml25_prepare_private_func(ifxc_mgga_func_type *p, ifxc_nspin nspin)
{
  memset(p, 0, sizeof(*p));
  p->info = &IFXC_ML25_PRIVATE_INFO;
  p->nspin = (nspin == IFXC_POLARIZED) ? IFXC_MGGA_POLARIZED : IFXC_MGGA_UNPOLARIZED;
  p->dens_threshold = p->info->dens_threshold;
  p->sigma_threshold = pow(p->dens_threshold, 4.0 / 3.0);
  p->zeta_threshold = DBL_EPSILON;
  p->tau_threshold = 1e-20;

  return IFXC_OK;
}

static void
ifxc_ml25_finish_private_func(ifxc_mgga_func_type *p)
{
  (void)p;
}

ifxc_status
ifxc_ml25_eval(
    const ifxc_func_type *func,
    const ifxc_input *input,
    size_t nentries,
    const ifxc_deriv_entry *entries)
{
  const ifxc_handle_impl *impl = (const ifxc_handle_impl *)((func != NULL) ? func->impl : NULL);

  if(impl == NULL || input == NULL || entries == NULL){
    return IFXC_E_INVALID_ARGUMENT;
  }
  if(impl->feature_set != IFXC_FEATURE_SET_ML25){
    return IFXC_E_UNKNOWN_FEATURE_SET;
  }
  if(nentries == 1 && entries[0].target == IFXC_TARGET_INTEGRAL &&
     entries[0].order == 0){
    return ifxc_ml25_eval_integral_order0(impl, input, entries[0].out);
  }
  {
    unsigned int max_order = ifxc_vector_mgga_max_requested_order(nentries, entries);
    if(max_order <= 3 && ifxc_ml25_has_combined_order(max_order)){
      return ifxc_ml25_eval_combined(impl, input, nentries, entries, max_order);
    }
    return IFXC_E_UNSUPPORTED_DERIVATIVE;
  }
}

static ifxc_status
ifxc_ml25_eval_combined(
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
  size_t npoints;
  ifxc_mgga_funcs selected_work;
  double feature_scale = 1.0;

  npoints = input->npoints;
  ifxc_vector_mgga_set_private_dimensions((impl->nspin == IFXC_UNPOLARIZED) ? IFXC_MGGA_UNPOLARIZED : IFXC_MGGA_POLARIZED,
                                          &private_dims);
  ifxc_vector_mgga_zero_dimensions_above_order(&private_dims, max_order);
  ifxc_vector_mgga_scale_output_dimensions(&private_dims, impl->nfeatures);

  status = ifxc_vector_mgga_transpose_inputs(&impl->dims, input, &rho_tm, &sigma_tm, &lapl_tm, &tau_tm);
  if(status != IFXC_OK){
    goto cleanup;
  }

  status = ifxc_vector_mgga_allocate_private_outputs(&private_dims, npoints, &out);
  if(status != IFXC_OK){
    goto cleanup;
  }

  status = ifxc_ml25_prepare_private_func(&p, impl->nspin);
  if(status != IFXC_OK){
    goto cleanup;
  }
  p.dim = private_dims;

  selected_work = (impl->nspin == IFXC_UNPOLARIZED)
    ? ifxc_ml25_combined_work_mgga()->unpol[max_order]
    : ifxc_ml25_combined_work_mgga()->pol[max_order];
  if(selected_work == NULL){
    status = IFXC_E_INTERNAL;
    goto cleanup;
  }

  selected_work(&p, npoints, rho_tm, sigma_tm, lapl_tm, tau_tm, &out);

  for(entry_index = 0; entry_index < nentries; ++entry_index){
    const ifxc_deriv_entry *entry = &entries[entry_index];
    const double *src = ifxc_vector_mgga_select_source_buffer(&out, entry);

    if(src == NULL){
      status = IFXC_E_INTERNAL;
      goto cleanup;
    }
    status = ifxc_vector_mgga_copy_combined_entry(
        input, &impl->dims, impl->nfeatures, entry, src, feature_scale);
    if(status != IFXC_OK){
      goto cleanup;
    }
  }

  status = IFXC_OK;

cleanup:
  ifxc_ml25_finish_private_func(&p);
  ifxc_vector_mgga_free_private_outputs(&out);
  ifxc_vector_mgga_free_inputs(rho_tm, sigma_tm, lapl_tm, tau_tm);
  return status;
}

ifxc_status
ifxc_ml25_eval_first_derivatives_contracted(
    const ifxc_handle_impl *impl,
    const ifxc_input *input,
    const double *coeffs,
    double *d_rho,
    double *d_sigma,
    double *d_tau)
{
  ifxc_mgga_dimensions private_dims;
  ifxc_mgga_func_type p = {0};
  ifxc_status status;

  if(impl == NULL || input == NULL || coeffs == NULL ||
     d_rho == NULL || d_sigma == NULL || d_tau == NULL){
    return IFXC_E_INVALID_ARGUMENT;
  }
  if(impl->feature_set != IFXC_FEATURE_SET_ML25){
    return IFXC_E_UNKNOWN_FEATURE_SET;
  }
  if(!ifxc_ml25_has_combined_order(1)){
    return IFXC_E_UNSUPPORTED_DERIVATIVE;
  }

  ifxc_vector_mgga_set_private_dimensions(
      (impl->nspin == IFXC_UNPOLARIZED) ? IFXC_MGGA_UNPOLARIZED : IFXC_MGGA_POLARIZED,
      &private_dims);
  ifxc_vector_mgga_zero_dimensions_above_order(&private_dims, 1);
  private_dims.zk = 0;

  status = ifxc_ml25_prepare_private_func(&p, impl->nspin);
  if(status != IFXC_OK){
    goto cleanup;
  }
  p.dim = private_dims;
  p.feature_coeffs = coeffs;

  status = ifxc_ml25_combined_eval_contracted_component_major(
      &p, input->npoints, input->rho, input->sigma, input->tau,
      d_rho, d_sigma, d_tau);

cleanup:
  ifxc_ml25_finish_private_func(&p);
  return status;
}

ifxc_status
ifxc_ml25_eval_integral_order0(
    const ifxc_handle_impl *impl,
    const ifxc_input *input,
    double *features)
{
  ifxc_mgga_dimensions private_dims;
  ifxc_mgga_func_type p = {0};
  ifxc_status status;

  if(impl == NULL || input == NULL || features == NULL ||
     input->weights == NULL){
    return IFXC_E_INVALID_ARGUMENT;
  }
  if(impl->feature_set != IFXC_FEATURE_SET_ML25){
    return IFXC_E_UNKNOWN_FEATURE_SET;
  }

  ifxc_vector_mgga_set_private_dimensions(
      (impl->nspin == IFXC_UNPOLARIZED)
          ? IFXC_MGGA_UNPOLARIZED
          : IFXC_MGGA_POLARIZED,
      &private_dims);
  ifxc_vector_mgga_zero_dimensions_above_order(&private_dims, 0);
  ifxc_vector_mgga_scale_output_dimensions(&private_dims, impl->nfeatures);
  status = ifxc_ml25_prepare_private_func(&p, impl->nspin);
  if(status != IFXC_OK){
    return status;
  }
  p.dim = private_dims;
  status = ifxc_ml25_combined_integral_order0_component_major(
      &p, input->npoints, input->rho, input->sigma, input->tau,
      input->weights, features);
  ifxc_ml25_finish_private_func(&p);
  return status;
}
