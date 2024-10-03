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

// void HDF5_IO_Write_0_opensbliblock00(ops_block& opensbliblock00, ops_dat& phi_B0, ops_dat& x0_B0, int HDF5_timing){
// double cpu_start0, elapsed_start0;
// if (HDF5_timing == 1){
// ops_timers(&cpu_start0, &elapsed_start0);
// }
// // Writing OPS datasets
// char name0[80];
// sprintf(name0, "opensbli_output.h5");
// ops_fetch_block_hdf5_file(opensbliblock00, name0);
// ops_fetch_dat_hdf5_file(phi_B0, name0);
// ops_fetch_dat_hdf5_file(x0_B0, name0);
// // Writing simulation constants
// write_constants(name0);
// if (HDF5_timing == 1){
// double cpu_end0, elapsed_end0;
// ops_timers(&cpu_end0, &elapsed_end0);
// ops_printf("-----------------------------------------\n");
// ops_printf("Time to write HDF5 file: %s: %lf\n", name0, elapsed_end0-elapsed_start0);
// ops_printf("-----------------------------------------\n");

// }
// }

#define NAME_LENGTH 80

ops_dat take_strided_slice(ops_dat dat, int stride) {

  for (int i = 0; i < OPS_MAX_DIM; ++i) {
    dat->size[i] = ceil(dat->size[i] / stride);
  }

  /* these actually need to be the sum of size, rather than indexed to [0] */
  int new_size[OPS_MAX_DIM] = {dat->size[0] - dat->d_p[0] + dat->d_m[0]};
  char *data_slice = (char *)malloc(new_size[0] * sizeof(char));

  /* So the problem here is that ops_dat->data is a char*, so we need to be
     smarter about how we copy the data from ops_data->data to data and how
     we set the type */
  /* We also should start from the halo and not from the beginning, and vice
     versa for the end of the array */
  for (int i = abs(dat->d_m[0]); i < new_size[0] - dat->d_p[0]; i++) { /* again, needs to over the sum of new_size */
    const int stride_index = i * stride;
    if (stride_index >= block0np0) {
      break;
    }
    data_slice[i] = dat->data[stride_index];
  }

  // return ops_decl_dat(dat->block, 1, new_size, dat->base, dat->d_m, dat->d_p, data_slice, dat->type, dat->name);
  return ops_decl_dat(dat->block, dat->dim, dat->size, dat->base, dat->d_m, dat->d_p, data_slice, dat->type,
                      "phi_B0_strided");
}

void hdf5_strided(ops_block &block, ops_dat &phi_B0, ops_dat &x0_B0, int stride, int hdf5_timing) {
  char name[NAME_LENGTH] = "opensbli_output.h5";

  double cpu_start0, elapsed_start0;
  if (hdf5_timing == 1) {
    ops_timers(&cpu_start0, &elapsed_start0);
  }

  // Taking strided slice of each ops_dat dataset
  ops_dat phi_B0_strided = take_strided_slice(phi_B0, stride);

  // Writing OPS datasets
  ops_fetch_block_hdf5_file(block, name);
  // ops_fetch_dat_hdf5_file(phi_B0, name);
  ops_fetch_dat_hdf5_file(phi_B0_strided, name);
  ops_fetch_dat_hdf5_file(x0_B0, name);

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
