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
 * Takes an input dataset to create a new strided dataset
 */
ops_dat take_strided_slice(ops_dat dat, int stride) {
  /*
   * First, figure out the size of the dataset minus the padding for halos and
   * memory alignment (x_pad, I think). We need to re-create d_p without padding
   * because it contains the value of xpad in the original ops_dat.
   */
  int d_p[OPS_MAX_DIM] = {0};
  int new_size[OPS_MAX_DIM] = {1};
  for (int i = 0; i < OPS_MAX_DIM; ++i) {
    if (dat->size[i] > 1) {
      /*
       * Could instead use block0np0 / stride, because we know during code
       * generation the size of the block
       */
      new_size[i] = (dat->size[i] + dat->d_m[i] - dat->d_p[i]) / stride;
    }
  }
  for (int i = 0; i < OPS_MAX_DIM; ++i) {
    d_p[i] = dat->d_p[i] - dat->x_pad;
  }
  /*
   * Create a new data set, where double *_dummy is required so ops_decl_dat
   * will allocate enough memory as it takes the type of dummy to determine the
   * size of each element in the data set.
   */
  double *_dummy = NULL; /* With code generation, this is needs to be smarter */
  ops_dat new_dat =
      ops_decl_dat(dat->block, dat->dim, new_size, dat->base, dat->d_m, d_p, _dummy, dat->type, dat->name);
  /*
   * Using memcpy, copy the ops_dat->data into the new ops_dat->data, one
   * element at a time. The indices are important, as we want to ignore the
   * halos, I think.
   */
  size_t position = 0;
  for (int i = 0; i < OPS_MAX_DIM; ++i) {
    for (int j = 0; j < dat->size[i]; j += stride) {
      const size_t index_old = j * dat->elem_size;
      const size_t index_new = position * dat->elem_size;
      memcpy(&new_dat->data[index_new], &dat->data[index_old], dat->elem_size);
      position += 1;
    }
  }

  return new_dat;
}

/*
 * Main steering function for strided HDF5 output
 */
void hdf5_strided(ops_block &block, ops_dat &phi_B0, ops_dat &x0_B0, int stride, int hdf5_timing) {
  char name[80] = "opensbli_output-strided.h5";

  double cpu_start0, elapsed_start0;
  if (hdf5_timing == 1) {
    ops_timers(&cpu_start0, &elapsed_start0);
  }

  // Taking strided slice of each ops_dat dataset
  ops_dat phi_B0_strided = take_strided_slice(phi_B0, stride);
  ops_dat x0_B0_strided = take_strided_slice(x0_B0, stride);

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
