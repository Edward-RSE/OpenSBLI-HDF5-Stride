/*
 * Writes constants to the HDF5
 */
void write_constants(const char *filename) {
  ops_write_const_hdf5("Delta0block0", 1, "double", (char *)&Delta0block0, filename);
  ops_write_const_hdf5("HDF5_timing", 1, "int", (char *)&HDF5_timing, filename);
  ops_write_const_hdf5("block0np0", 1, "int", (char *)&block0np0, filename);
  ops_write_const_hdf5("c0", 1, "double", (char *)&c0, filename);
  ops_write_const_hdf5("dt", 1, "double", (char *)&dt, filename);
  ops_write_const_hdf5("niter", 1, "int", (char *)&niter, filename);
  ops_write_const_hdf5("simulation_time", 1, "double", (char *)&simulation_time, filename);
  ops_write_const_hdf5("start_iter", 1, "int", (char *)&start_iter, filename);
  ops_write_const_hdf5("iter", 1, "int", (char *)&iter, filename);
}

/*
 * Original steering function for outputting HDF5 data
 */
void HDF5_IO_Write_0_opensbliblock00(ops_block &opensbliblock00, ops_dat &phi_B0, ops_dat &x0_B0, int HDF5_timing) {
  double cpu_start0, elapsed_start0;
  if (HDF5_timing == 1) {
    ops_timers(&cpu_start0, &elapsed_start0);
  }
  // Writing OPS datasets
  char name0[80] = "opensbli_output.h5";
  ops_fetch_block_hdf5_file(opensbliblock00, name0);
  ops_fetch_dat_hdf5_file(phi_B0, name0);
  ops_fetch_dat_hdf5_file(x0_B0, name0);
  // Writing simulation constants
  write_constants(name0);
  if (HDF5_timing == 1) {
    double cpu_end0, elapsed_end0;
    ops_timers(&cpu_end0, &elapsed_end0);
    ops_printf("-----------------------------------------\n");
    ops_printf("Time to write HDF5 file: %s: %lf\n", name0, elapsed_end0 - elapsed_start0);
    ops_printf("-----------------------------------------\n");
  }
}

/*
 * Main steering function for strided HDF5 output
 */

#include "io_strided.h"

void HDF5_IO_Write_Strided_0_opensbliblock00(ops_block &block, ops_dat &phi_B0, ops_dat &x0_B0, int hdf5_timing) {
  char name[80] = "opensbli_output-strided.h5";

  double cpu_start0, elapsed_start0;
  if (hdf5_timing == 1) {
    ops_timers(&cpu_start0, &elapsed_start0);
  }

  // Taking strided slice of each ops_dat dataset
  const size_t stride_i = 10;
  ops_dat x0_B0_strided = create_strided_ops_dat_with_halo_cells(x0_B0, stride_i);
  ops_dat phi_B0_strided = create_strided_ops_dat_with_halo_cells(phi_B0, stride_i);

  // Writing OPS datasets
  ops_fetch_block_hdf5_file(block, name);
  ops_fetch_dat_hdf5_file(x0_B0_strided, name);
  ops_fetch_dat_hdf5_file(phi_B0_strided, name);

  // Writing simulation constants
  write_constants(name);

  if (hdf5_timing == 1) {
    double cpu_end0, elapsed_end0;
    ops_timers(&cpu_end0, &elapsed_end0);
    ops_printf("-----------------------------------------\n");
    ops_printf("Time to write HDF5 file: %s: %lf\n", name, elapsed_end0 - elapsed_start0);
    ops_printf("-----------------------------------------\n");
  }
}
