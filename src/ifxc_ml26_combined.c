#include "util.h"
#include "ifxc.h"

#include <string.h>

#define XC_PREFIX ml26_combined_
#define XC_MAPLE2C_INC "maple2c/if_mgga/mgga_xc_ml26.c"
#include "mgga_xc_ml25_inc.c"

ifxc_mgga_funcs_variants *
ifxc_ml26_combined_work_mgga(void)
{
  return &ml26_combined_work_mgga;
}

int
ifxc_ml26_has_combined_order(unsigned int order)
{
  if(order >= 5){
    return 0;
  }
  return ml26_combined_work_mgga.unpol[order] != NULL &&
         ml26_combined_work_mgga.pol[order] != NULL;
}

static int
ifxc_ml26_prepare_component_major_point(
    const ifxc_mgga_func_type *p,
    size_t npoints,
    size_t point,
    const double *rho,
    const double *sigma,
    const double *tau,
    double local_rho[2],
    double local_sigma[3],
    double local_tau[2])
{
  double density = rho[point];

  if(p->nspin == IFXC_MGGA_POLARIZED){
    density += rho[npoints + point];
  }
  if(density < p->dens_threshold){
    return 0;
  }

  local_rho[0] = m_max(p->dens_threshold, rho[point]);
  local_sigma[0] = m_max(
      p->sigma_threshold * p->sigma_threshold, sigma[point]);
  local_tau[0] = m_max(p->tau_threshold, tau[point]);

  if(p->nspin == IFXC_MGGA_POLARIZED){
    double sigma_average;
    local_rho[1] = m_max(p->dens_threshold, rho[npoints + point]);
    local_sigma[2] = m_max(
        p->sigma_threshold * p->sigma_threshold,
        sigma[2 * npoints + point]);
    local_tau[1] = m_max(p->tau_threshold, tau[npoints + point]);
    sigma_average = 0.5 * (local_sigma[0] + local_sigma[2]);
    local_sigma[1] = sigma[npoints + point];
    local_sigma[1] = m_max(-sigma_average, local_sigma[1]);
    local_sigma[1] = m_min(+sigma_average, local_sigma[1]);
  }

  return 1;
}

ifxc_status
ifxc_ml26_combined_eval_contracted_component_major(
    const ifxc_mgga_func_type *p,
    size_t npoints,
    const double *rho,
    const double *sigma,
    const double *tau,
    double *d_rho,
    double *d_sigma,
    double *d_tau)
{
  size_t point;

  if(p == NULL || d_rho == NULL || d_sigma == NULL || d_tau == NULL ||
     p->feature_coeffs == NULL){
    return IFXC_E_INVALID_ARGUMENT;
  }
  if(npoints == 0){
    return IFXC_OK;
  }
  if(rho == NULL || sigma == NULL || tau == NULL){
    return IFXC_E_INVALID_ARGUMENT;
  }

  for(point = 0; point < npoints; ++point){
    double local_rho[2] = {0.0, 0.0};
    double local_sigma[3] = {0.0, 0.0, 0.0};
    double local_tau[2] = {0.0, 0.0};
    double local_lapl[2] = {0.0, 0.0};
    double local_d_rho[2] = {0.0, 0.0};
    double local_d_sigma[3] = {0.0, 0.0, 0.0};
    double local_d_tau[2] = {0.0, 0.0};
    ifxc_mgga_out_params out = {0};
    size_t component;

    if(!ifxc_ml26_prepare_component_major_point(
           p, npoints, point, rho, sigma, tau,
           local_rho, local_sigma, local_tau)){
      for(component = 0; component < (size_t)p->dim.vrho; ++component){
        d_rho[component * npoints + point] = 0.0;
      }
      for(component = 0; component < (size_t)p->dim.vsigma; ++component){
        d_sigma[component * npoints + point] = 0.0;
      }
      for(component = 0; component < (size_t)p->dim.vtau; ++component){
        d_tau[component * npoints + point] = 0.0;
      }
      continue;
    }

    out.vrho = local_d_rho;
    out.vsigma = local_d_sigma;
    out.vtau = local_d_tau;
    if(p->nspin == IFXC_MGGA_UNPOLARIZED){
      ml26_combined_func_vxc_unpol(
          p, 0, local_rho, local_sigma, local_lapl, local_tau, &out);
    }else{
      ml26_combined_func_vxc_pol(
          p, 0, local_rho, local_sigma, local_lapl, local_tau, &out);
    }

    for(component = 0; component < (size_t)p->dim.vrho; ++component){
      d_rho[component * npoints + point] = local_d_rho[component];
    }
    for(component = 0; component < (size_t)p->dim.vsigma; ++component){
      d_sigma[component * npoints + point] = local_d_sigma[component];
    }
    for(component = 0; component < (size_t)p->dim.vtau; ++component){
      d_tau[component * npoints + point] = local_d_tau[component];
    }
  }

  return IFXC_OK;
}

ifxc_status
ifxc_ml26_combined_eval_order0_component_major(
    const ifxc_mgga_func_type *p,
    size_t npoints,
    const double *rho,
    const double *sigma,
    const double *tau,
    const double *weights,
    double *local_values,
    double *features)
{
  size_t point;

  if(p == NULL || (local_values == NULL && features == NULL)){
    return IFXC_E_INVALID_ARGUMENT;
  }
  if(features != NULL){
    if(weights == NULL){
      return IFXC_E_INVALID_ARGUMENT;
    }
    memset(features, 0, IFXC_ML26_NFEATURES * sizeof(double));
  }
  if(npoints == 0){
    return IFXC_OK;
  }
  if(rho == NULL || sigma == NULL || tau == NULL){
    return IFXC_E_INVALID_ARGUMENT;
  }

  for(point = 0; point < npoints; ++point){
    double local_rho[2] = {0.0, 0.0};
    double local_sigma[3] = {0.0, 0.0, 0.0};
    double local_tau[2] = {0.0, 0.0};
    double local_lapl[2] = {0.0, 0.0};
    double generated[IFXC_ML26_NFEATURES] = {0.0};
    ifxc_mgga_out_params out = {0};
    double density;
    size_t feature;

    if(!ifxc_ml26_prepare_component_major_point(
           p, npoints, point, rho, sigma, tau,
           local_rho, local_sigma, local_tau)){
      if(local_values != NULL){
        memset(local_values + point * IFXC_ML26_NFEATURES, 0,
               IFXC_ML26_NFEATURES * sizeof(double));
      }
      continue;
    }

    out.zk = generated;
    if(p->nspin == IFXC_MGGA_UNPOLARIZED){
      ml26_combined_func_exc_unpol(
          p, 0, local_rho, local_sigma, local_lapl, local_tau, &out);
      density = rho[point];
    }else{
      ml26_combined_func_exc_pol(
          p, 0, local_rho, local_sigma, local_lapl, local_tau, &out);
      density = rho[point] + rho[npoints + point];
    }

    for(feature = 0; feature < IFXC_ML26_NFEATURES; ++feature){
      double value = density * generated[feature];
      if(local_values != NULL){
        local_values[point * IFXC_ML26_NFEATURES + feature] = value;
      }
      if(features != NULL){
        features[feature] += weights[point] * value;
      }
    }
  }

  return IFXC_OK;
}
