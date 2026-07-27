#ifndef IFXC_H
#define IFXC_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define IFXC_VERSION_MAJOR 0
#define IFXC_VERSION_MINOR 4
#define IFXC_VERSION_PATCH 0
#define IFXC_API_VERSION 4

#define IFXC_ML25_NFEATURES 66
#define IFXC_ML26_NFEATURES 69

#define IFXC_ML25_LAK_X 0
#define IFXC_ML25_LAK_C 1
#define IFXC_ML25_LYP_T1 2
#define IFXC_ML25_LYP_T2 3
#define IFXC_ML25_LYP_T3 4
#define IFXC_ML25_LYP_T4 5
#define IFXC_ML25_LYP_T5 6
#define IFXC_ML25_LYP_T6 7
#define IFXC_ML25_MN15_C01 8
#define IFXC_ML25_MN15_C02 9
#define IFXC_ML25_MN15_C03 10
#define IFXC_ML25_MN15_C04 11
#define IFXC_ML25_MN15_C05 12
#define IFXC_ML25_MN15_C06 13
#define IFXC_ML25_MN15_C07 14
#define IFXC_ML25_MN15_C08 15
#define IFXC_ML25_MN15_C09 16
#define IFXC_ML25_MN15_C10 17
#define IFXC_ML25_MN15_C11 18
#define IFXC_ML25_MN15_C12 19
#define IFXC_ML25_MN15_C13 20
#define IFXC_ML25_MN15_C14 21
#define IFXC_ML25_MN15_C15 22
#define IFXC_ML25_MN15_C16 23
#define IFXC_ML25_MN15_C17 24
#define IFXC_ML25_MN15_C18 25
#define IFXC_ML25_MN15_C19 26
#define IFXC_ML25_MN15_C20 27
#define IFXC_ML25_MN15_C21 28
#define IFXC_ML25_MN15_C22 29
#define IFXC_ML25_MN15_C23 30
#define IFXC_ML25_MN15_C24 31
#define IFXC_ML25_MN15_C25 32
#define IFXC_ML25_MN15_C26 33
#define IFXC_ML25_MN15_C27 34
#define IFXC_ML25_MN15_C28 35
#define IFXC_ML25_MN15_C29 36
#define IFXC_ML25_MN15_C30 37
#define IFXC_ML25_MN15_C31 38
#define IFXC_ML25_MN15_C32 39
#define IFXC_ML25_MN15_C33 40
#define IFXC_ML25_MN15_C34 41
#define IFXC_ML25_MN15_C35 42
#define IFXC_ML25_MN15_C36 43
#define IFXC_ML25_MN15_C37 44
#define IFXC_ML25_MN15_C38 45
#define IFXC_ML25_MN15_C39 46
#define IFXC_ML25_MN15_C40 47
#define IFXC_ML25_MN15_A01 48
#define IFXC_ML25_MN15_A02 49
#define IFXC_ML25_MN15_A03 50
#define IFXC_ML25_MN15_A04 51
#define IFXC_ML25_MN15_A05 52
#define IFXC_ML25_MN15_A06 53
#define IFXC_ML25_MN15_A07 54
#define IFXC_ML25_MN15_A08 55
#define IFXC_ML25_MN15_A09 56
#define IFXC_ML25_MN15_B01 57
#define IFXC_ML25_MN15_B02 58
#define IFXC_ML25_MN15_B03 59
#define IFXC_ML25_MN15_B04 60
#define IFXC_ML25_MN15_B05 61
#define IFXC_ML25_MN15_B06 62
#define IFXC_ML25_MN15_B07 63
#define IFXC_ML25_MN15_B08 64
#define IFXC_ML25_MN15_B09 65

/* ML26 preserves the complete ML25 vector at indices 0 through 65. */
#define IFXC_ML26_ML25_NFEATURES IFXC_ML25_NFEATURES
#define IFXC_ML26_CS1_SAME_SPIN_1 66
#define IFXC_ML26_CS1_SAME_SPIN_2 67
#define IFXC_ML26_CS1_OPPOSITE_SPIN 68

#if defined(_WIN32) || defined(__CYGWIN__)
#  if defined(IFXC_BUILD_DLL)
#    define IFXC_API __declspec(dllexport)
#  else
#    define IFXC_API __declspec(dllimport)
#  endif
#elif defined(__GNUC__) || defined(__clang__)
#  define IFXC_API __attribute__((visibility("default")))
#else
#  define IFXC_API
#endif

typedef enum {
  IFXC_OK = 0,
  IFXC_E_INVALID_ARGUMENT = 1,
  IFXC_E_ALLOCATION = 2,
  IFXC_E_UNKNOWN_FEATURE_SET = 3,
  IFXC_E_INVALID_NSPIN = 4,
  IFXC_E_MISSING_INPUT = 5,
  IFXC_E_UNSUPPORTED_DERIVATIVE = 6,
  IFXC_E_UNSUPPORTED_VARIABLE = 7,
  IFXC_E_NOT_INITIALIZED = 8,
  IFXC_E_INTERNAL = 9
} ifxc_status;

typedef enum {
  IFXC_UNPOLARIZED = 1,
  IFXC_POLARIZED = 2
} ifxc_nspin;

typedef enum {
  IFXC_FEATURE_SET_ML25 = 1,
  IFXC_FEATURE_SET_ML26 = 2
} ifxc_feature_set;

typedef enum {
  IFXC_VAR_RHO = 0,
  IFXC_VAR_SIGMA = 1,
  IFXC_VAR_LAPL = 2,
  IFXC_VAR_TAU = 3
} ifxc_variable;

typedef enum {
  IFXC_TARGET_LOCAL = 0,
  IFXC_TARGET_INTEGRAL = 1
} ifxc_target;

typedef enum {
  IFXC_FEATURE_KIND_EXCHANGE = 1,
  IFXC_FEATURE_KIND_CORRELATION = 2,
  IFXC_FEATURE_KIND_EXCHANGE_CORRELATION = 3
} ifxc_feature_kind;

/* The metadata queries use the plain names, so the struct aliases carry a
 * _t suffix to keep the C identifier namespace valid. */
struct ifxc_feature_set_info {
  int feature_set;
  const char *key;
  const char *name;
  size_t nfeatures;
  unsigned int max_deriv_order;
};

typedef struct ifxc_feature_set_info ifxc_feature_set_info_t;

struct ifxc_feature_info {
  int index;
  int paper_id;
  const char *key;
  const char *group;
  ifxc_feature_kind kind;
};

typedef struct ifxc_feature_info ifxc_feature_info_t;

struct ifxc_dimensions {
  size_t rho;
  size_t sigma;
  size_t lapl;
  size_t tau;
};

typedef struct ifxc_dimensions ifxc_dimensions_t;

typedef struct ifxc_func_type {
  void *impl;
} ifxc_func_type;

/* Input arrays are component-major:
 *   rho[component * npoints + point]
 *   sigma[component * npoints + point]
 *   tau[component * npoints + point]
 * The current feature sets do not use lapl; pass NULL for lapl. weights is
 * required only for entries whose target is IFXC_TARGET_INTEGRAL. */
typedef struct {
  size_t npoints;
  const double *rho;
  const double *sigma;
  const double *lapl;
  const double *tau;
  const double *weights;
} ifxc_input;

typedef struct {
  ifxc_target target;
  unsigned int order;
  /* Length is order. For order 0 this may be NULL. For order > 0, variables
   * must be sorted in nondecreasing enum order. */
  const ifxc_variable *vars;
  double *out;
} ifxc_deriv_entry;

/* Output layout:
 *   LOCAL, order 0:
 *     out[point * nfeatures + feature] = h_f(r_point)
 *   LOCAL, order > 0:
 *     out[(point * nfeatures + feature) * ncomponents + component]
 *       = d^order h_f(r_point) / dvars...
 *   INTEGRAL, order 0:
 *     out[feature] = sum_point weights[point] * h_f(r_point)
 *   INTEGRAL, order > 0:
 *     out[(point * nfeatures + feature) * ncomponents + component]
 *       = weights[point] * d^order h_f(r_point) / dvars...
 *
 * Repeated variables use compact symmetric component enumeration; mixed
 * variable groups use products of those compact groups. Use ifxc_output_size()
 * to size buffers for a requested entry; ifxc_output_size() validates entry
 * metadata but does not require entry->out to be set. ifxc_eval() requires all
 * requested entries to have non-NULL output buffers. */

IFXC_API const char *ifxc_version_string(void);
IFXC_API int ifxc_version_major(void);
IFXC_API int ifxc_version_minor(void);
IFXC_API int ifxc_version_patch(void);

IFXC_API const char *ifxc_strerror(ifxc_status status);

IFXC_API ifxc_status ifxc_feature_set_info(
    int feature_set,
    const ifxc_feature_set_info_t **info);

IFXC_API ifxc_status ifxc_feature_info(
    int feature_set,
    size_t index,
    const ifxc_feature_info_t **info);

IFXC_API ifxc_status ifxc_dimensions(
    int feature_set,
    ifxc_nspin nspin,
    ifxc_dimensions_t *dims);

IFXC_API ifxc_status ifxc_init(
    ifxc_func_type *func,
    int feature_set,
    ifxc_nspin nspin);

IFXC_API void ifxc_end(ifxc_func_type *func);

IFXC_API ifxc_status ifxc_nfeatures(
    const ifxc_func_type *func,
    size_t *nfeatures);

IFXC_API ifxc_status ifxc_max_deriv_order(
    const ifxc_func_type *func,
    unsigned int *max_order);

IFXC_API ifxc_status ifxc_func_dimensions(
    const ifxc_func_type *func,
    ifxc_dimensions_t *dims);

IFXC_API ifxc_status ifxc_output_size(
    const ifxc_func_type *func,
    const ifxc_input *input,
    const ifxc_deriv_entry *entry,
    size_t *n_double);

IFXC_API ifxc_status ifxc_eval(
    const ifxc_func_type *func,
    const ifxc_input *input,
    size_t nentries,
    const ifxc_deriv_entry *entries);

/* Fused first-derivative contraction for the feature set owned by func.
 *
 * coeffs has ifxc_nfeatures(func) entries. Outputs are component-major:
 *   d_rho[component * npoints + point]
 *   d_sigma[component * npoints + point]
 *   d_tau[component * npoints + point]
 *
 * This is equivalent to evaluating local first derivatives for rho, sigma,
 * and tau with ifxc_eval() and contracting each feature derivative with coeffs.
 * weights are not used. */
IFXC_API ifxc_status ifxc_eval_first_derivatives_contracted(
    const ifxc_func_type *func,
    const ifxc_input *input,
    const double *coeffs,
    double *d_rho,
    double *d_sigma,
    double *d_tau);

/* ABI-compatible ML25 name retained for existing callers. */
IFXC_API ifxc_status ifxc_eval_ml25_first_derivatives_contracted(
    const ifxc_func_type *func,
    const ifxc_input *input,
    const double *coeffs,
    double *d_rho,
    double *d_sigma,
    double *d_tau);

#ifdef __cplusplus
}
#endif

#endif /* IFXC_H */
