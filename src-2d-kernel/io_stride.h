
#include <ops_seq.h>

/**
 * @brief
 *
 * @param original_dat
 * @param strided_data
 */
void kernel_copy_to_strided_dat(const ACC<double> &original_dat, ACC<double> &strided_data) {
  strided_data(0, 0) = original_dat(0, 0);
}

/**
 * @brief
 *
 */

ops_dat rho_B0_strided;
ops_stencil stencil2d_00;
ops_stencil stencil2d_00_strided;

int stride[] = {2, 2};

void declare_strided_datasets(ops_block block) {
  int strided_size[] = {block0np0 / stride[0], block0np1 / stride[1]};
  int strided_base[] = {0, 0};
  int strided_d_p[] = {0, 0};
  int strided_d_m[] = {0, 0};
  double *dummy = NULL;

  rho_B0_strided =
      ops_decl_dat(block, 1, strided_size, strided_base, strided_d_m, strided_d_p, dummy, "double", "rho_B0_strided");

  int stencil_point[] = {0, 0};
  // int restrict_stencil_point[] = {0, 0, -1, 0, 1, 0};
  stencil2d_00 = ops_decl_stencil(2, 1, stencil_point, "stencil2d_00");
  stencil2d_00_strided = ops_decl_strided_stencil(2, 1, stencil_point, stride, "stencil2d_00_strided");
}

void write_strided_data_sets(ops_block block, ops_dat rho_B0) {
  char name[80] = "opensbli_output-strided.h5";
  int iter_range[] = {0, block0np0 / stride[0], 0, block0np1 / stride[1]};

  ops_printf("\nIter range %d %d %d %d\n", iter_range[0], iter_range[1], iter_range[2], iter_range[3]);

  ops_par_loop(kernel_copy_to_strided_dat, "kernel_copy_to_strided_dat", block, 2, iter_range,
               ops_arg_dat(rho_B0, 1, stencil2d_00_strided, "double", OPS_READ),
               ops_arg_dat(rho_B0_strided, 1, stencil2d_00, "double", OPS_WRITE));

  ops_fetch_block_hdf5_file(block, name);
  ops_fetch_dat_hdf5_file(rho_B0_strided, name);
}
