/**
 * @brief Function for writing strided data to disk
 */

#include <ops_seq.h> /* Have to do this for some reason... it's not been resolved properly */

/*
 * Declare datasets in global scope because I'm being lazy
 */
ops_dat rho_B0_strided;

/*
 * Declare stencils in global scope because I'm being lazy
 */
ops_stencil stencil2d_00;
ops_stencil stencil2d_00_strided;

/**
 * @brief Copy data from one dataset to another.
 *
 * @param original_dat  Source dataset
 * @param strided_data  Destination dataset
 */
void restrict_strided_copy_kernel(const ACC<double> &original_dat, ACC<double> &strided_data) {
  strided_data(0, 0) = original_dat(0, 0);
}

/**
 * @brief Control function for copying datasets
 *
 * @param block
 * @param stride
 * @param original_dat
 * @param strided_dat
 *
 * @details
 *
 * In this function, we use a parallel loop per dataset because that's how
 * the kernel is set up. However it would be possible (and maybe more efficient,
 * I haven't benchmarked it) to do all copies in a single parallel loop/kernel
 * because they will all have the same interation range and the same stencils.
 *
 */
void copy_to_strided_dat(ops_block block, int stride[], ops_dat &original_dat, ops_dat &strided_dat) {
  const int dims = 2;
  int iter_range[] = {0, block0np0 / stride[0], 0, block0np1 / stride[1]};
  /*
   * Use a parallel loop to copy data from the original to the smaller data
   * set. The important thing here is that we are looping over the smaller
   * range, e.g. block0np0 / stride[0] rather than block0np0 AND that we are
   * using a strided/restricted stencil for the larger dataset.
   */
  ops_par_loop(restrict_strided_copy_kernel, "restrict_strided_copy_kernel", block, dims, iter_range,
               ops_arg_dat(original_dat, 1, stencil2d_00_strided, "double", OPS_READ),
               ops_arg_dat(strided_dat, 1, stencil2d_00, "double", OPS_WRITE));
}

/**
 * @brief Initialise stencil and datasets for strided data output
 *
 * @param block   The OPS block datasets will be associated with
 * @param stride  The strided of the output datasets
 *
 * @details
 *
 * Note that THIS HAS TO BE CALLED BEFORE OPS_PARTITION, otherwise in MPI modes
 * the program will crash because it will not be able to partition the datasets
 * across ranks.
 *
 */
void HDF5_IO_Init_0_opensbliblock00_strided(ops_block block, int stride[]) {
  /*
   * Determine the size of the strided datasets. We have some halo cells, but we
   * don't actually need them if we don't want them.
   */
  int strided_size[] = {block0np0 / stride[0], block0np1 / stride[1]};
  int strided_base[] = {0, 0};
  int strided_d_p[] = {5, 5};
  int strided_d_m[] = {-5, -5};
  double *dummy = NULL;

  /*
   * Declare smaller datasets to store strided versions.
   */
  rho_B0_strided =
      ops_decl_dat(block, 1, strided_size, strided_base, strided_d_m, strided_d_p, dummy, "double", "rho_B0_strided");

  /*
   * Declare TWO stencils. The first stencil is a regular 2D stencil whilst the
   * second is a strided/restricted 2D stencil. Both stencils are around the
   * point (0, 0) as we don't require any other adjacent cells to copy data.
   */
  const int dims = 2;
  const int points = 1;
  int stencil_point[] = {0, 0};
  stencil2d_00 = ops_decl_stencil(dims, points, stencil_point, "stencil2d_00");
  stencil2d_00_strided = ops_decl_restrict_stencil(dims, points, stencil_point, stride, "RESTRICT_stencil2d_00");
}

/**
 * @brief Write data in a strided output to disk
 *
 * @param block   The OPS block containing the datasets
 * @param stride  The strided of the output datasets
 * @param rho_B0  The first dataset to write
 */
void HDF5_IO_Write_0_opensbliblock00_strided(ops_block block, int stride[], ops_dat &rho_B0) {
  double cpu_start0;
  double cpu_end0;
  double elapsed_start0;
  double elapsed_end0;
  const char name[80] = "opensbli_output-strided.h5";

  ops_timers(&cpu_start0, &elapsed_start0);

  /* Copy data to strided datasets */
  copy_to_strided_dat(block, stride, rho_B0, rho_B0_strided);

  /* Write to disk */
  ops_fetch_block_hdf5_file(block, name);
  ops_fetch_dat_hdf5_file(rho_B0_strided, name);

  ops_timers(&cpu_end0, &elapsed_end0);
  ops_printf("-----------------------------------------\n");
  ops_printf("Time to write strided HDF5 file: %s: %lf\n", name, elapsed_end0 - elapsed_start0);
  ops_printf("-----------------------------------------\n");
}
