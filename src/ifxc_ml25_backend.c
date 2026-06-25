#include "ifxc_internal.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "xc.h"

extern xc_mgga_funcs_variants *work_mgga_ml25[];

static const xc_func_info_type IFXC_ML25_PRIVATE_INFO = {
  .number = 888,
  .flags = XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | XC_FLAGS_HAVE_ALL,
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

static void ifxc_ml25_free_inputs(double *rho_tm, double *sigma_tm, double *lapl_tm, double *tau_tm);
static void ifxc_ml25_free_private_outputs(xc_mgga_out_params *out);

static void
ifxc_ml25_set_private_dimensions(int nspin, xc_dimensions *dim)
{
  memset(dim, 0, sizeof(*dim));
  dim->zk = 1;

  if(nspin == XC_UNPOLARIZED){
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

static const double *
ifxc_ml25_select_source_buffer(const xc_mgga_out_params *out, const ifxc_deriv_entry *entry)
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

static ifxc_status
ifxc_ml25_copy_feature_entry(
    const ifxc_input *input,
    const ifxc_dimensions_t *dims,
    size_t feature_index,
    size_t nfeatures,
    const ifxc_deriv_entry *entry,
    const double *src,
    double scale)
{
  size_t npoints = input->npoints;
  size_t ncomp = 0;
  size_t point;
  ifxc_status status = ifxc_derivative_component_count(dims, entry, &ncomp);
  if(status != IFXC_OK){
    return status;
  }

  if(entry->target == IFXC_TARGET_LOCAL){
    if(entry->order == 0){
      for(point = 0; point < npoints; ++point){
        entry->out[point * nfeatures + feature_index] = scale * src[point];
      }
      return IFXC_OK;
    }

    {
      size_t comp;
      for(point = 0; point < npoints; ++point){
        for(comp = 0; comp < ncomp; ++comp){
          size_t dest_index = (point * nfeatures + feature_index) * ncomp + comp;
          entry->out[dest_index] = scale * src[point * ncomp + comp];
        }
      }
    }
    return IFXC_OK;
  }

  if(entry->target == IFXC_TARGET_INTEGRAL){
    if(entry->order == 0){
      double sum = 0.0;
      for(point = 0; point < npoints; ++point){
        sum += input->weights[point] * scale * src[point];
      }
      entry->out[feature_index] = sum;
      return IFXC_OK;
    }

    {
      size_t comp;
      for(point = 0; point < npoints; ++point){
        for(comp = 0; comp < ncomp; ++comp){
          size_t dest_index = (point * nfeatures + feature_index) * ncomp + comp;
          entry->out[dest_index] = input->weights[point] * scale * src[point * ncomp + comp];
        }
      }
    }
    return IFXC_OK;
  }

  return IFXC_E_INVALID_ARGUMENT;
}

static ifxc_status
ifxc_ml25_transpose_inputs(
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
    ifxc_ml25_free_inputs(*rho_tm, *sigma_tm, *lapl_tm, *tau_tm);
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

static void
ifxc_ml25_free_inputs(double *rho_tm, double *sigma_tm, double *lapl_tm, double *tau_tm)
{
  free(rho_tm);
  free(sigma_tm);
  free(lapl_tm);
  free(tau_tm);
}

static ifxc_status
ifxc_ml25_allocate_private_outputs(
    const xc_dimensions *dims,
    size_t npoints,
    xc_mgga_out_params *out)
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
  ifxc_ml25_free_private_outputs(out);
  memset(out, 0, sizeof(*out));
  return status;
}

static void
ifxc_ml25_free_private_outputs(xc_mgga_out_params *out)
{
#define IFXC_FREE_FIELD(name) free(out->name);

  IFXC_ML25_PRIVATE_FIELDS(IFXC_FREE_FIELD);

#undef IFXC_FREE_FIELD
}

static void
ifxc_ml25_zero_private_outputs(
    const xc_dimensions *dims,
    size_t npoints,
    xc_mgga_out_params *out)
{
#define IFXC_ZERO_FIELD(name) \
  do { \
    if(out->name != NULL && npoints > 0 && dims->name > 0){ \
      memset(out->name, 0, npoints * (size_t)dims->name * sizeof(double)); \
    } \
  } while(0);

  IFXC_ML25_PRIVATE_FIELDS(IFXC_ZERO_FIELD)

#undef IFXC_ZERO_FIELD
}

static ifxc_status
ifxc_ml25_prepare_private_func(xc_func_type *p, ifxc_nspin nspin)
{
  memset(p, 0, sizeof(*p));
  p->info = (xc_func_info_type *)&IFXC_ML25_PRIVATE_INFO;
  p->nspin = (nspin == IFXC_POLARIZED) ? XC_POLARIZED : XC_UNPOLARIZED;
  p->dens_threshold = p->info->dens_threshold;
  p->sigma_threshold = pow(p->dens_threshold, 4.0 / 3.0);
  p->zeta_threshold = DBL_EPSILON;
  p->tau_threshold = 1e-20;

  return IFXC_OK;
}

static void
ifxc_ml25_finish_private_func(xc_func_type *p)
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
  xc_dimensions private_dims;
  double *rho_tm = NULL;
  double *sigma_tm = NULL;
  double *lapl_tm = NULL;
  double *tau_tm = NULL;
  xc_mgga_out_params out = {0};
  xc_func_type p = {0};
  ifxc_status status;
  size_t feature_index;
  size_t entry_index;
  size_t npoints;
  size_t eval_npoints;
  xc_mgga_funcs selected_work;
  double feature_scale;

  if(impl == NULL || input == NULL || entries == NULL){
    return IFXC_E_INVALID_ARGUMENT;
  }
  if(impl->feature_set != IFXC_FEATURE_SET_ML25){
    return IFXC_E_UNKNOWN_FEATURE_SET;
  }

  npoints = input->npoints;
  feature_scale = 2.0;
  eval_npoints = npoints;
  ifxc_ml25_set_private_dimensions((impl->nspin == IFXC_UNPOLARIZED) ? XC_UNPOLARIZED : XC_POLARIZED,
                                   &private_dims);

  status = ifxc_ml25_transpose_inputs(&impl->dims, input, &rho_tm, &sigma_tm, &lapl_tm, &tau_tm);
  if(status != IFXC_OK){
    goto cleanup;
  }

  status = ifxc_ml25_allocate_private_outputs(&private_dims, npoints, &out);
  if(status != IFXC_OK){
    goto cleanup;
  }

  status = ifxc_ml25_prepare_private_func(&p, impl->nspin);
  if(status != IFXC_OK){
    goto cleanup;
  }
  p.dim = private_dims;

  for(feature_index = 0; feature_index < impl->nfeatures; ++feature_index){
    ifxc_ml25_zero_private_outputs(&private_dims, npoints, &out);
    selected_work = (impl->nspin == IFXC_UNPOLARIZED)
      ? work_mgga_ml25[feature_index + 1]->unpol[4]
      : work_mgga_ml25[feature_index + 1]->pol[4];
    if(selected_work == NULL){
      status = IFXC_E_INTERNAL;
      goto cleanup;
    }
    selected_work(&p, eval_npoints, rho_tm, sigma_tm, lapl_tm, tau_tm, &out);

    for(entry_index = 0; entry_index < nentries; ++entry_index){
      const ifxc_deriv_entry *entry = &entries[entry_index];
      const double *src = ifxc_ml25_select_source_buffer(&out, entry);

      if(src == NULL){
        status = IFXC_E_INTERNAL;
        goto cleanup;
      }
      status = ifxc_ml25_copy_feature_entry(
          input, &impl->dims, feature_index, impl->nfeatures, entry, src, feature_scale);
      if(status != IFXC_OK){
        goto cleanup;
      }
    }
  }

  status = IFXC_OK;

cleanup:
  ifxc_ml25_finish_private_func(&p);
  ifxc_ml25_free_private_outputs(&out);
  ifxc_ml25_free_inputs(rho_tm, sigma_tm, lapl_tm, tau_tm);
  return status;
}
