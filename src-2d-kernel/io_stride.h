/**
 * @brief Function for writing strided data to disk
 */

/*
 * Declare datasets in global scope because I'm being lazy
 */
ops_dat rho_B0_strided;

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

  ops_printf("Strided size: %d %d\n", strided_size[0], strided_size[1]);
  ops_printf("Stride size: %d %d\n", stride[0], stride[1]);

  /*
   * Declare smaller datasets to store strided versions.
   */
  rho_B0_strided = ops_decl_dat(block, 1, strided_size, strided_base, strided_d_m, strided_d_p, stride, dummy, "double",
                                "rho_B0_strided");

  /*
   * Declare TWO stencils. The first stencil is a regular 2D stencil whilst the
   * second is a strided/restricted 2D stencil. Both stencils are around the
   * point (0, 0) as we don't require any other adjacent cells to copy data.
   */
  const int dims = 2;
  const int points = 1;
  int s2D_00[] = {0, 0};
  int s2D_00_M10_P10[] = {0, 0, -1, 0, 1, 0};
  S2D_00 = ops_decl_stencil(2, 1, s2D_00, "00");
  S2D_RESTRICT_00_M10_P10 = ops_decl_restrict_stencil(2, 3, s2D_00_M10_P10, stride, "RESTRICT_00_M10_P10");
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
