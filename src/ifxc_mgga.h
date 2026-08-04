/* SPDX-License-Identifier: MPL-2.0 */

#ifndef IFXC_MGGA_H
#define IFXC_MGGA_H

#include <stddef.h>

#define IFXC_MGGA_UNPOLARIZED 1
#define IFXC_MGGA_POLARIZED 2

#define IFXC_MGGA_FLAGS_HAVE_EXC  (1 << 0)
#define IFXC_MGGA_FLAGS_HAVE_VXC  (1 << 1)
#define IFXC_MGGA_FLAGS_HAVE_FXC  (1 << 2)
#define IFXC_MGGA_FLAGS_HAVE_KXC  (1 << 3)
#define IFXC_MGGA_FLAGS_HAVE_LXC  (1 << 4)
#define IFXC_MGGA_FLAGS_3D        (1 << 7)
#define IFXC_MGGA_FLAGS_NEEDS_LAPLACIAN (1 << 15)
#define IFXC_MGGA_FLAGS_NEEDS_TAU (1 << 16)
#define IFXC_MGGA_FLAGS_ENFORCE_FHC (1 << 17)

#define IFXC_MGGA_FLAGS_HAVE_ALL \
  (IFXC_MGGA_FLAGS_HAVE_EXC | IFXC_MGGA_FLAGS_HAVE_VXC | \
   IFXC_MGGA_FLAGS_HAVE_FXC | IFXC_MGGA_FLAGS_HAVE_KXC | \
   IFXC_MGGA_FLAGS_HAVE_LXC)

typedef struct ifxc_mgga_func_info {
  int flags;
  double dens_threshold;
} ifxc_mgga_func_info;

typedef struct ifxc_mgga_dimensions {
  int rho;
  int sigma;
  int lapl;
  int tau;

  int zk;
  int vrho;
  int vsigma;
  int vlapl;
  int vtau;

  int v2rho2;
  int v2rhosigma;
  int v2rholapl;
  int v2rhotau;
  int v2sigma2;
  int v2sigmalapl;
  int v2sigmatau;
  int v2lapl2;
  int v2lapltau;
  int v2tau2;

  int v3rho3;
  int v3rho2sigma;
  int v3rho2lapl;
  int v3rho2tau;
  int v3rhosigma2;
  int v3rhosigmalapl;
  int v3rhosigmatau;
  int v3rholapl2;
  int v3rholapltau;
  int v3rhotau2;
  int v3sigma3;
  int v3sigma2lapl;
  int v3sigma2tau;
  int v3sigmalapl2;
  int v3sigmalapltau;
  int v3sigmatau2;
  int v3lapl3;
  int v3lapl2tau;
  int v3lapltau2;
  int v3tau3;

  int v4rho4;
  int v4rho3sigma;
  int v4rho3lapl;
  int v4rho3tau;
  int v4rho2sigma2;
  int v4rho2sigmalapl;
  int v4rho2sigmatau;
  int v4rho2lapl2;
  int v4rho2lapltau;
  int v4rho2tau2;
  int v4rhosigma3;
  int v4rhosigma2lapl;
  int v4rhosigma2tau;
  int v4rhosigmalapl2;
  int v4rhosigmalapltau;
  int v4rhosigmatau2;
  int v4rholapl3;
  int v4rholapl2tau;
  int v4rholapltau2;
  int v4rhotau3;
  int v4sigma4;
  int v4sigma3lapl;
  int v4sigma3tau;
  int v4sigma2lapl2;
  int v4sigma2lapltau;
  int v4sigma2tau2;
  int v4sigmalapl3;
  int v4sigmalapl2tau;
  int v4sigmalapltau2;
  int v4sigmatau3;
  int v4lapl4;
  int v4lapl3tau;
  int v4lapl2tau2;
  int v4lapltau3;
  int v4tau4;
} ifxc_mgga_dimensions;

typedef struct ifxc_mgga_out_params {
  double *zk;
  double *vrho;
  double *vsigma;
  double *vlapl;
  double *vtau;
  double *v2rho2;
  double *v2rhosigma;
  double *v2rholapl;
  double *v2rhotau;
  double *v2sigma2;
  double *v2sigmalapl;
  double *v2sigmatau;
  double *v2lapl2;
  double *v2lapltau;
  double *v2tau2;
  double *v3rho3;
  double *v3rho2sigma;
  double *v3rho2lapl;
  double *v3rho2tau;
  double *v3rhosigma2;
  double *v3rhosigmalapl;
  double *v3rhosigmatau;
  double *v3rholapl2;
  double *v3rholapltau;
  double *v3rhotau2;
  double *v3sigma3;
  double *v3sigma2lapl;
  double *v3sigma2tau;
  double *v3sigmalapl2;
  double *v3sigmalapltau;
  double *v3sigmatau2;
  double *v3lapl3;
  double *v3lapl2tau;
  double *v3lapltau2;
  double *v3tau3;
  double *v4rho4;
  double *v4rho3sigma;
  double *v4rho3lapl;
  double *v4rho3tau;
  double *v4rho2sigma2;
  double *v4rho2sigmalapl;
  double *v4rho2sigmatau;
  double *v4rho2lapl2;
  double *v4rho2lapltau;
  double *v4rho2tau2;
  double *v4rhosigma3;
  double *v4rhosigma2lapl;
  double *v4rhosigma2tau;
  double *v4rhosigmalapl2;
  double *v4rhosigmalapltau;
  double *v4rhosigmatau2;
  double *v4rholapl3;
  double *v4rholapl2tau;
  double *v4rholapltau2;
  double *v4rhotau3;
  double *v4sigma4;
  double *v4sigma3lapl;
  double *v4sigma3tau;
  double *v4sigma2lapl2;
  double *v4sigma2lapltau;
  double *v4sigma2tau2;
  double *v4sigmalapl3;
  double *v4sigmalapl2tau;
  double *v4sigmalapltau2;
  double *v4sigmatau3;
  double *v4lapl4;
  double *v4lapl3tau;
  double *v4lapl2tau2;
  double *v4lapltau3;
  double *v4tau4;
} ifxc_mgga_out_params;

struct ifxc_mgga_func_type;

typedef void (*ifxc_mgga_funcs)(
    const struct ifxc_mgga_func_type *p,
    size_t np,
    const double *rho,
    const double *sigma,
    const double *lapl,
    const double *tau,
    ifxc_mgga_out_params *out);

typedef struct ifxc_mgga_funcs_variants {
  const ifxc_mgga_funcs unpol[5];
  const ifxc_mgga_funcs pol[5];
} ifxc_mgga_funcs_variants;

typedef struct ifxc_mgga_func_type {
  const ifxc_mgga_func_info *info;
  int nspin;
  ifxc_mgga_dimensions dim;
  const double *feature_coeffs;
  double dens_threshold;
  double zeta_threshold;
  double sigma_threshold;
  double tau_threshold;
} ifxc_mgga_func_type;

#endif /* IFXC_MGGA_H */
