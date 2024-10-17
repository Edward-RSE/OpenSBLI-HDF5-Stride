void write_constants(const char *filename) {
  ops_write_const_hdf5("DRP_filt", 1, "double", (char *)&DRP_filt, filename);
  ops_write_const_hdf5("Delta0block0", 1, "double", (char *)&Delta0block0, filename);
  ops_write_const_hdf5("Delta1block0", 1, "double", (char *)&Delta1block0, filename);
  ops_write_const_hdf5("Delta2block0", 1, "double", (char *)&Delta2block0, filename);
  ops_write_const_hdf5("HDF5_timing", 1, "int", (char *)&HDF5_timing, filename);
  ops_write_const_hdf5("Minf", 1, "double", (char *)&Minf, filename);
  ops_write_const_hdf5("Pr", 1, "double", (char *)&Pr, filename);
  ops_write_const_hdf5("Re", 1, "double", (char *)&Re, filename);
  ops_write_const_hdf5("block0np0", 1, "int", (char *)&block0np0, filename);
  ops_write_const_hdf5("block0np1", 1, "int", (char *)&block0np1, filename);
  ops_write_const_hdf5("block0np2", 1, "int", (char *)&block0np2, filename);
  ops_write_const_hdf5("dt", 1, "double", (char *)&dt, filename);
  ops_write_const_hdf5("filter_frequency", 1, "int", (char *)&filter_frequency, filename);
  ops_write_const_hdf5("gama", 1, "double", (char *)&gama, filename);
  ops_write_const_hdf5("niter", 1, "int", (char *)&niter, filename);
  ops_write_const_hdf5("simulation_time", 1, "double", (char *)&simulation_time, filename);
  ops_write_const_hdf5("start_iter", 1, "int", (char *)&start_iter, filename);
  ops_write_const_hdf5("write_output_file", 1, "int", (char *)&write_output_file, filename);
  ops_write_const_hdf5("iter", 1, "int", (char *)&iter, filename);
}

void HDF5_IO_Write_0_opensbliblock00_dynamic(ops_block &opensbliblock00, int iter, ops_dat &rho_B0, ops_dat &rhou0_B0,
                                             ops_dat &rhou1_B0, ops_dat &rhou2_B0, ops_dat &rhoE_B0, int HDF5_timing) {
  double cpu_start0, elapsed_start0;
  if (HDF5_timing == 1) {
    ops_timers(&cpu_start0, &elapsed_start0);
  }
  // Writing OPS datasets
  char name0[80];
  sprintf(name0, "opensbli_output_%06d.h5", iter + 1);
  ops_fetch_block_hdf5_file(opensbliblock00, name0);
  ops_fetch_dat_hdf5_file(rho_B0, name0);
  ops_fetch_dat_hdf5_file(rhou0_B0, name0);
  ops_fetch_dat_hdf5_file(rhou1_B0, name0);
  ops_fetch_dat_hdf5_file(rhou2_B0, name0);
  ops_fetch_dat_hdf5_file(rhoE_B0, name0);
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

void HDF5_IO_Write_0_opensbliblock00(ops_block &opensbliblock00, ops_dat &rho_B0, ops_dat &rhou0_B0, ops_dat &rhou1_B0,
                                     ops_dat &rhou2_B0, ops_dat &rhoE_B0, int HDF5_timing) {
  double cpu_start0, elapsed_start0;
  if (HDF5_timing == 1) {
    ops_timers(&cpu_start0, &elapsed_start0);
  }
  // Writing OPS datasets
  char name0[80];
  sprintf(name0, "opensbli_output.h5");
  ops_fetch_block_hdf5_file(opensbliblock00, name0);
  ops_fetch_dat_hdf5_file(rho_B0, name0);
  ops_fetch_dat_hdf5_file(rhou0_B0, name0);
  ops_fetch_dat_hdf5_file(rhou1_B0, name0);
  ops_fetch_dat_hdf5_file(rhou2_B0, name0);
  ops_fetch_dat_hdf5_file(rhoE_B0, name0);
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
 * Takes an input dataset to create a new strided dataset
 * NOTE: this has been written in a style that it has been code generated
 */
ops_dat create_strided_ops_dat(ops_dat dat, size_t stride) {

  int d_p[] = {5, 5, 5};
  int d_m[] = {-5, -5, -5};
  int size[] = {block0np0 / stride, block0np1 / stride, block0np2 / stride};

  /*
   * Create a new data set, where double *_dummy is required so ops_decl_dat
   * will allocate enough memory as it takes the type of _dummy to determine the
   * size of each element in the data set.
   */
  double *_dummy = NULL;
  ops_dat new_dat = ops_decl_dat(dat->block, dat->dim, size, dat->base, d_m, d_p, _dummy, dat->type, dat->name);

  /*
   * `position` is the index where we are writing data to in the new dataset.
   * It has to start at the first element in the new dataset, ignoring the first
   * set of halo cells.
   */
  // size_t position = abs(dat->d_m[0]);
  size_t position = 0;

  /*
   * Copies one element at a time from one char* to another char*, hence why
   * there is a `memcpy` rather than straight up assignment. `position` is used
   * to track the index into the strided dataset (`new_dat`). Note that each
   * index variable is incremented by stride.
   *
   * We need to be careful about the indices for each loop, as we don't want to
   * copy the halo cells. This means we need to start each loop at abs(d_m[i])
   * and then add that amount to the other end.
   */

  int offset[] = {abs(d_m[0]), abs(d_m[1]), abs(d_m[2])};

  for (int k = offset[2]; k < size[2] + offset[2]; k += stride) {
    //
    for (int j = offset[1]; j < size[1] + offset[1]; j += stride) {
      //
      for (int i = offset[0]; i < size[0] + offset[0]; i += stride) {
        //
        const size_t index_dat = i + j * size[0] + k * size[0] * size[1];
        const size_t index_new = position * new_dat->elem_size;
        memcpy(new_dat->data + index_new, dat->data + index_dat, new_dat->elem_size);
        position += 1;

        printf("i: %d, j: %d, k: %d, position: %d, index_dat: %d, index_new: %d\n", i, j, k, position, index_dat,
               index_new);
      }
      position += offset[0] + d_p[0];
    }
    position += 2 * offset[1] * new_dat->size[0];
  }

  /*
   * Increment `position` so we don't write to space reserved for halos and
   * memory padding
   */
  // if (i < (size_t)OPS_MAX_DIM - 1) {
  //   position += abs(dat->d_m[i + 1]) + dat->d_p[i];
  // }

  return new_dat;
}

void HDF5_IO_Write_Strided_0_opensbliblock00(ops_block &opensbliblock00, ops_dat &rho_B0, ops_dat &rhou0_B0,
                                             ops_dat &rhou1_B0, ops_dat &rhou2_B0, ops_dat &rhoE_B0, int HDF5_timing) {
  double cpu_start0, elapsed_start0;
  if (HDF5_timing == 1) {
    ops_timers(&cpu_start0, &elapsed_start0);
  }

  // Writing OPS datasets
  char name0[80];
  sprintf(name0, "opensbli_output-strided.h5");

  int stride = 1;
  ops_dat rho_B0_strided = create_strided_ops_dat(rho_B0, stride);
  ops_dat rhou0_B0_strided = create_strided_ops_dat(rhou0_B0, stride);
  ops_dat rhou1_B0_strided = create_strided_ops_dat(rhou1_B0, stride);
  ops_dat rhou2_B0_strided = create_strided_ops_dat(rhou2_B0, stride);
  ops_dat rhoE_B0_strided = create_strided_ops_dat(rhoE_B0, stride);

  ops_fetch_block_hdf5_file(opensbliblock00, name0);
  ops_fetch_dat_hdf5_file(rho_B0_strided, name0);
  ops_fetch_dat_hdf5_file(rhou0_B0_strided, name0);
  ops_fetch_dat_hdf5_file(rhou1_B0_strided, name0);
  ops_fetch_dat_hdf5_file(rhou2_B0_strided, name0);
  ops_fetch_dat_hdf5_file(rhoE_B0_strided, name0);

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
