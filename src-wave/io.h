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
ops_dat take_strided_slice(ops_dat dat, size_t stride) {
  /*
   * First, figure out the size of the dataset minus the padding for halos and
   * memory alignment (x_pad, I think). We need to re-create d_p without padding
   * because it contains the value of xpad in the original ops_dat.
   */
  int d_p[dat->dim];
  int new_size[dat->dim];
  for (size_t i = 0; i < (size_t)dat->dim; ++i) {
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
   * will allocate enough memory as it takes the type of _dummy to determine the
   * size of each element in the data set.
   */
  double *_dummy = NULL; /* With code generation, this is needs to be smarter */
  ops_dat new_dat =
      ops_decl_dat(dat->block, dat->dim, new_size, dat->base, dat->d_m, d_p, _dummy, dat->type, dat->name);

  /*
   * `position` is the index where we are writing data to in the new dataset.
   * It has to start at the first element in the new dataset, ignoring the first
   * set of halo cells.
   */
  size_t position = abs(dat->d_m[0]);
  for (size_t i = 0; i < (size_t)dat->dim; ++i) {
    /*
     * Start is the first element after the "negative" halos and stop is, or
     * should, be the last element before the "positive" halos
     */
    const size_t start = abs(dat->d_m[i]);
    const size_t stop = dat->size[i] + start;

    /*
     * Copies one element at a time from one char* to another char*, hence why
     * there is a memcpy rather than straight up assignment. `position` is used
     * to track the index into the strided dataset. Note that the index variable
     * is incremented by stride.
     */
    for (size_t j = start; j < stop; j += stride) {
      const size_t index_dat = j * dat->elem_size;
      const size_t index_new = position * dat->elem_size;
      memcpy(new_dat->data + index_new, dat->data + index_dat, dat->elem_size);
      position++;
    }
    /*
     * Increment `position` so we don't write to space reserved for halos and
     * memory padding
     */
    if (i < (size_t)dat->dim - 1) {
      position += abs(dat->d_m[i + 1]) + dat->d_p[i];
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
