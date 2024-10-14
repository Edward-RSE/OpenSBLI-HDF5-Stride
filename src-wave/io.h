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
   * First, figure out he size of the dataset minus the padding for halos
   */
  int total_size = 0;
  int new_size[OPS_MAX_DIM] = {1};
  for (int i = 0; i < OPS_MAX_DIM; ++i) {
    total_size += dat->size[i];
    new_size[i] = dat->size[i]; //- dat->d_p[i] + dat->d_m[i];
  }

  /* Next, we need to get data from dat->data into a new array. We have to
   * essentially do this element-by-element and be careful of padding. This
   * padding includes MPI halos and padding for memory alignment
   */
  size_t new_data_position = 0;
  char *new_data = (char *)ops_malloc(total_size * dat->elem_size);

  /*
   * Let's get this working in 1D for now, and disregard anything to do with
   * other dimensions. So this effectively removes the loop over OPS_MAX_DIM and
   * instead we'll loop over the dat->size[0] excluding the padding
   */

  char name[80];
  snprintf(name, 80, "_debug/%s-%d.txt", dat->name, stride);
  FILE *fp = fopen(name, "w");

  double sneak_peak = -1;
  for (int i = 0; i < dat->size[0]; ++i) {
    const size_t index_new = new_data_position * dat->elem_size;
    const size_t index_dat = i * dat->elem_size;
    memcpy(&new_data[index_new], &dat->data[index_dat], dat->elem_size);

    /*
     * This is for debug
     */
    memcpy(&sneak_peak, &new_data[index_new], dat->elem_size);
    fprintf(fp, "%e\n", sneak_peak);

    new_data_position += 1;
  }

  fclose(fp);

  /*
   * Create a new data set and then replace the old data and size with the
   * new ones we've created
   */
  double *_dummy = NULL;
  ops_dat new_dat =
      ops_decl_dat(dat->block, dat->dim, dat->size, dat->base, dat->d_m, dat->d_p, _dummy, dat->type, dat->name);
  for (int i = 0; i < OPS_MAX_DIM; ++i) {
    new_dat->size[i] = new_size[i];
  }
  new_dat->data = new_data;
  new_dat->mem = total_size * new_dat->elem_size;

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
