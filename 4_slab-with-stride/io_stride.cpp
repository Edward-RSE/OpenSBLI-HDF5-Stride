/**
 * @file io_stride.cpp
 * @author Edward Parkinson
 * @brief Implementation of strided IO
 * @date November 2024
 *
 */

#define OPS_3D
#define OPS_API 2

#include "ops_seq.h"

/*
 * Declare datasets and stencils in global scope because I'm being lazy
 */
ops_dat rho_B0_strided;
ops_stencil stencil2d_00;
ops_stencil stencil2d_restrict_00;

/**
 * @brief Copy data from one dataset to another.
 *
 * @param original_dat  Source dataset
 * @param strided_dat   Destination dataset
 */
void restrict_kernel(const ACC<double> &original_dat, ACC<double> &strided_dat, const int *idx) {
  strided_dat(0, 0, 0) = original_dat(0, 0, 0);
}

/**
 * @brief Control function for copying datasets
 *
 * @param block         The OPS block datasets are associated with
 * @param block0np0     The size of the original dataset in direction 0
 * @param block0np1     The size of the original dataset in direction 1
 * @param stride        The stride of the output strided dataset
 * @param original_dat  The original dataset
 * @param strided_dat   The (smaller) strided dataset to copy data to
 *
 * @details
 *
 * In this function, we use a parallel loop per dataset because that's how
 * the kernel is set up. However it would be possible (and maybe more efficient,
 * I haven't benchmarked it) to do all copies in a single parallel loop/kernel
 * because they will all have the same interation range and the same stencils.
 *
 */
void copy_to_strided_dat(ops_block block, int iter_range[], ops_dat &original_dat, ops_dat &strided_dat) {
  /*
   * Use a parallel loop to copy data from the original to the smaller data
   * set. The important thing here is that we are looping over the smaller
   * range, e.g. block0np0 / stride[0] rather than block0np0 AND that we are
   * using a strided/restricted stencil for the larger dataset.
   */
  ops_par_loop(restrict_kernel, "restrict_kernel", block, 3, iter_range,
               ops_arg_dat(original_dat, 1, stencil2d_restrict_00, "double", OPS_READ),
               ops_arg_dat(strided_dat, 1, stencil2d_00, "double", OPS_WRITE), ops_arg_idx());
}

/**
 * @brief Initialise stencil and datasets for strided data output
 *
 * @param block      The OPS block datasets will be associated with
 * @param block0np0  The size of the original dataset in direction 0
 * @param block0np1  The size of the original dataset in direction 1
 * @param stride     The strided of the output datasets
 *
 * @details
 *
 * Note that THIS HAS TO BE CALLED BEFORE OPS_PARTITION, otherwise in MPI modes
 * the program will crash because it will not be able to partition the datasets
 * across ranks.
 *
 */
void HDF5_IO_Init_0_opensbliblock00_strided(ops_block block, int block0np0, int block0np1, int block0np2,
                                            int stride[]) {
  /*
   * Declare the slab range and size, following logic described in the comment
   * around line 130 in `HDF5_IO_Write_0_opensbliblock00_slab`.
   */
  const int offset0 = 40;
  const int offset1 = 40;
  const int offset2 = 10;
  int slab_size[] = {2 * offset0, 2 * offset1, 2 * offset2};

  /*
   * Determine the size of the strided datasets. We have some halo cells, but we
   * don't actually need them if we don't want them.
   */
  int strided_size[] = {block0np0 / stride[0], block0np1 / stride[1], block0np2 / stride[2]};
  int strided_base[] = {0, 0, 0};
  int strided_d_p[] = {5, 5, 5};
  int strided_d_m[] = {-5, -5, -5};
  double *dummy = NULL;

  /*
   * Declare smaller datasets to store strided versions.
   */
  rho_B0_strided = ops_decl_dat(block, 1, strided_size, strided_base, strided_d_m, strided_d_p, stride, dummy, "double",
                                "rho_B0_strided");

  /*
   * Declare TWO stencils. The first stencil is a regular 2D stencil whilst the
   * second is a strided/restricted 2D stencil. Both stencils are around the
   * point (0, 0) as we don't require any other adjacent cells to copy data.
   * The restrict stencil has to have the word RESTRICT in it, otherwise this
   * will not work.
   */
  int stencil_point[] = {0, 0, 0};
  stencil2d_00 = ops_decl_stencil(3, 1, stencil_point, "stencil_00");
  stencil2d_restrict_00 = ops_decl_restrict_stencil(3, 1, stencil_point, stride, "stencil_RESTRICT_00");
}

/**
 * @brief Write data in a strided output to disk
 *
 * @param block   The OPS block containing the datasets
 * @param stride  The strided of the output datasets
 * @param rho_B0  The first dataset to write
 */
void HDF5_IO_Write_0_opensbliblock00_slab(char name[], ops_block block, int block0np0, int block0np1, int block0np2,
                                          int stride[], ops_dat &rho_B0) {
  double cpu_start0;
  double cpu_end0;
  double elapsed_start0;
  double elapsed_end0;
  char output_name[80] = "opensbli_output-slab_strided.h5";

  /* Declare slab range - offset here is referring to the number of cells on
   * either side of some central point we're taking the slab from. For example
   * imagine we have offset0 = 5. The method assumes (although it doesn't have
   * to) that the slab will start at some point X - offset and will span to
   * X + offset where X is some arbitrary point to define the slab.
   */
  const int offset0 = 40;
  const int offset1 = 40;
  const int offset2 = 10;
  const int X = block0np0 / 2;
  const int Y = block0np1 / 2;
  const int Z = block0np2 / 2;

  /* These are the coordinates of the slab we want to create */
  int slab_range[] = {X - offset0, X + offset0, Y - offset1, Y + offset1, Z - offset2, Z + offset2};

  /*
   * These are the coordinates of the strided slab we want to create - this needs
   * to be used in the parallel loop to copy data in conjuction with the strided
   * stencil
   */
  int stride_range[] = {slab_range[0], slab_range[1] / stride[0], slab_range[2], slab_range[3] / stride[1],
                        slab_range[4], slab_range[5] / stride[2]};

  /* Not sure if required */
  // int output_range[] = {0, 2 * offset0 / stride[0], 0, 2 * offset1 / stride[1], 0, 2 * offset2 / stride[2]};

  ops_timers(&cpu_start0, &elapsed_start0);

  /* Copy data to strided datasets */
  copy_to_strided_dat(block, stride_range, rho_B0, rho_B0_strided);

  /* Write to disk, using the standard HDF5 API */
  // char slab_name[80];
  // snprintf(slab_name, 80, "%s/rho_B0", name);
  // ops_write_data_slab_hdf5(rho_B0_strided, output_range, output_name, slab_name);
  ops_fetch_block_hdf5_file(block, output_name);
  ops_fetch_dat_hdf5_file(rho_B0_strided, output_name);

  ops_timers(&cpu_end0, &elapsed_end0);
  ops_printf("-----------------------------------------\n");
  ops_printf("Time to write HDF5 slab: %s: %lf\n", name, elapsed_end0 - elapsed_start0);
  ops_printf("-----------------------------------------\n");
}
