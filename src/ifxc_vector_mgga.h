/* SPDX-License-Identifier: MPL-2.0 */

#ifndef IFXC_VECTOR_MGGA_H
#define IFXC_VECTOR_MGGA_H

#include "ifxc_internal.h"
#include "ifxc_mgga.h"

void ifxc_vector_mgga_set_private_dimensions(
    int nspin,
    ifxc_mgga_dimensions *dim);

const double *ifxc_vector_mgga_select_source_buffer(
    const ifxc_mgga_out_params *out,
    const ifxc_deriv_entry *entry);

void ifxc_vector_mgga_zero_dimensions_above_order(
    ifxc_mgga_dimensions *dims,
    unsigned int max_order);

void ifxc_vector_mgga_scale_output_dimensions(
    ifxc_mgga_dimensions *dims,
    size_t nfeatures);

ifxc_status ifxc_vector_mgga_copy_combined_entry(
    const ifxc_input *input,
    const ifxc_dimensions_t *dims,
    size_t nfeatures,
    const ifxc_deriv_entry *entry,
    const double *src,
    double scale);

ifxc_status ifxc_vector_mgga_transpose_inputs(
    const ifxc_dimensions_t *dims,
    const ifxc_input *input,
    double **rho_tm,
    double **sigma_tm,
    double **lapl_tm,
    double **tau_tm);

void ifxc_vector_mgga_free_inputs(
    double *rho_tm,
    double *sigma_tm,
    double *lapl_tm,
    double *tau_tm);

ifxc_status ifxc_vector_mgga_allocate_private_outputs(
    const ifxc_mgga_dimensions *dims,
    size_t npoints,
    ifxc_mgga_out_params *out);

void ifxc_vector_mgga_free_private_outputs(ifxc_mgga_out_params *out);

unsigned int ifxc_vector_mgga_max_requested_order(
    size_t nentries,
    const ifxc_deriv_entry *entries);

#endif /* IFXC_VECTOR_MGGA_H */
