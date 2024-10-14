FILE *_debugfp;

void _debug_init(const char *name) {
  char fn[80];
  snprintf(fn, 80, "_debug/%s.txt", name);
  _debugfp = fopen(fn, "w");
}

void _debug_close(void) { fclose(_debugfp); }

void _debug_fprint(const char *element) {
  double value = -1;
  memcpy(&value, element, sizeof(double));
  fprintf(_debugfp, "%e\n", value);
}

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
  _debug_init(dat->name);
  /*
   * First, figure out the size of the dataset minus the padding for halos and
   * memory alignment (x_pad, I think). We also need to know how much elements
   * there are in total, so we can allocate enough memory for new_data.
   */
  size_t total_size = 0;
  int new_size[OPS_MAX_DIM] = {-1};
  for (int i = 0; i < OPS_MAX_DIM; ++i) {
    if (dat->size[i] > 1) {
      new_size[i] = (dat->size[i] + dat->d_m[i] - dat->d_p[i]) / stride;
    } else {
      new_size[i] = dat->size[i];
    }
    total_size += new_size[i];
  }

  /*
   * Create a new data set and then replace the old data and size with the
   * new ones we've created. double *_dummy is required so ops_decl_dat will
   * allocate enough memory as it takes the type of dummy to determine the size
   * of each element in the data set.
   */
  double *_dummy = NULL;
  ops_dat new_dat =
      ops_decl_dat(dat->block, dat->dim, new_size, dat->base, dat->d_m, dat->d_p, _dummy, dat->type, dat->name);

  size_t position = 0;

  for (int i = 0; i < dat->size[0]; i += stride) {
    const size_t index_old = i * dat->elem_size;
    const size_t index_new = position * dat->elem_size;
    memcpy(&new_dat->data[index_new], &dat->data[index_old], dat->elem_size);
    position += 1;
    _debug_fprint(&new_dat->data[index_new]);
  }

  _debug_close();

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
