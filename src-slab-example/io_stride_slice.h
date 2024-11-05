
/**
 * @brief Convert between 3D indices and linear indices
 *
 * @param[in] i The i index
 * @param[in] j The j index
 * @param[in] k The k index
 * @param[in] size The size of each dimension
 *
 * @returns The linear index
 *
 */
size_t index3d(size_t i, size_t j, size_t k, int *size) { return i + j * size[0] + k * size[0] * size[1]; }

/*
 * Declare global variables to use for strided IO
 */
ops_dat x0_B0_strided;
ops_dat x1_B0_strided;
ops_dat rho_B0_strided;
ops_dat rhou0_B0_strided;
ops_dat rhou1_B0_strided;
ops_dat rhou2_B0_strided;
ops_dat rhoE_B0_strided;
ops_dat WENO_filter_B0_strided;

void declare_empty_strided_ops_dat(ops_block block, int stride[2]) {
  /*
   * Determine the size of the strided dat halos and the actual size of the
   * dataset. In this case, we do not have any halo cells.
   */
  int strided_d_p[] = {0, 0, 0};
  int strided_d_m[] = {0, 0, 0};
  int strided_base[] = {0, 0, 0};
  int strided_size[] = {
      (int)ceil(block0np0 / stride[0]),
      (int)ceil(block0np1 / stride[1]),
      block0np2,
  };

  /*
   * Create a new dataset using the sizes above. We need to create a new ops_dat
   * rather than just a buffer since all the HDF5 functions (which we want to
   * re-use) expect an ops_dat and not a buffer containing data. Note that we
   * have a dummy variable. This variable *HAS* to be the type of data which
   * will be stored in the ops_dat. If it is not, this will allocate an
   * incorrect amount of data and will also initialise some of the fields, such
   * as ops_dat->elem_size, incorrectly. Anything wrong there will result in
   * copy errors and garbage output.
   */
  double *dummy1 = NULL, *dummy2 = NULL, *dummy3 = NULL, *dummy4 = NULL, *dummy5 = NULL, *dummy6 = NULL, *dummy7 = NULL,
         *dummy8 = NULL;

  x0_B0_strided =
      ops_decl_dat(block, 1, strided_size, strided_base, strided_d_m, strided_d_p, dummy1, "double", "x0_B0");
  x1_B0_strided =
      ops_decl_dat(block, 1, strided_size, strided_base, strided_d_m, strided_d_p, dummy2, "double", "x1_B0");
  rho_B0_strided =
      ops_decl_dat(block, 1, strided_size, strided_base, strided_d_m, strided_d_p, dummy3, "double", "rho_B0");
  rhou0_B0_strided =
      ops_decl_dat(block, 1, strided_size, strided_base, strided_d_m, strided_d_p, dummy4, "double", "rhou0_B0");
  rhou1_B0_strided =
      ops_decl_dat(block, 1, strided_size, strided_base, strided_d_m, strided_d_p, dummy5, "double", "rhou1_B0");
  rhou2_B0_strided =
      ops_decl_dat(block, 1, strided_size, strided_base, strided_d_m, strided_d_p, dummy6, "double", "rhou2_B0");
  rhoE_B0_strided =
      ops_decl_dat(block, 1, strided_size, strided_base, strided_d_m, strided_d_p, dummy7, "double", "rhoE_B0");
  WENO_filter_B0_strided =
      ops_decl_dat(block, 1, strided_size, strided_base, strided_d_m, strided_d_p, dummy8, "double", "WENO_filter_B0");
}

/**
 * @brief Creates a strided ops_dat, without halo cells, from an existing
 *        ops_dat.
 *
 * @param[in] dat The ops_dat to create a strided ops_dat from
 * @param[out] strided_dat The strided ops_dat with data copied to it
 * @param[in] stride The stride in each dimension
 *
 * @details
 *
 * This function creates a strided data set by copying data from a ops_dat to
 * a strided ops_dat (e.g. a ops_dat which is smaller).
 *
 */
void copy_to_strided_dat(ops_dat original_dat, ops_dat strided_dat, int stride[3]) {
  ops_get_data(original_dat);
  /*
   * Loop variables to control the start and stop of the loop over each
   * dimension. These values are based on the size of the original dataset. The
   * loops start at the first non-halo cell and stop at the last non-halo cell.
   * The loop ranges are over the entire dataset, original_dat->size.
   */
  int original_dat_start[] = {(int)abs(original_dat->d_m[0]), (int)abs(original_dat->d_m[1]),
                              (int)abs(original_dat->d_m[2])};
  int original_dat_stop[] = {original_dat->size[0] - original_dat->d_p[0], original_dat->size[1] - original_dat->d_p[1],
                             original_dat->size[2] - original_dat->d_p[2]};
  int offset_in_strided_buf[] = {(int)abs(strided_dat->d_m[0]), (int)abs(strided_dat->d_m[1]),
                                 (int)abs(strided_dat->d_m[2])};

  /*
   * Loop over each dimension, with each loop variable being incremented by the
   * stride for that dimension.
   */
  const int k = (original_dat->size[2] - original_dat->d_p[2] + original_dat->d_m[2]) / 2;

  for (int j = original_dat_start[1]; j < original_dat_stop[1]; j += stride[1]) {
    for (int i = original_dat_start[0]; i < original_dat_stop[0]; i += stride[0]) {
      /*
       * Using `index3d` we can calculate the index in the flat 1d buffers
       * given the 3d coordinates. The coordinate transformation for the
       * strided index is rather complex and ugly to look at. The index is
       * multiplied by the element size of ops_dat, because ops_dat->data is
       * a char*, so we are working with bytes rather than structured data.
       */
      const size_t stride_i = (i - original_dat_start[0]) / stride[0] + offset_in_strided_buf[0];
      const size_t stride_j = (j - original_dat_start[1]) / stride[0] + offset_in_strided_buf[1];
      const size_t stride_k = k;

      const size_t index_dat = index3d(i, j, k, original_dat->size) * original_dat->elem_size;
      const size_t index_stride = index3d(stride_i, stride_j, stride_k, strided_dat->size) * strided_dat->elem_size;

      /*
       * Copy the data from the original dataset into the strided buffer. This
       * has to be done using memcpy because ops_dat->data is a char* so we
       * are working in bytes, rather than structured data.
       */
      memcpy(&strided_dat->data[index_stride], &original_dat->data[index_dat], strided_dat->elem_size);
    }
  }
}

void ops_write_plane_group_strided_coords(char *slice_name0, int stride[], ops_dat x0_B0, ops_dat x1_B0) {
  copy_to_strided_dat(x0_B0, x0_B0_strided, stride);
  copy_to_strided_dat(x1_B0, x1_B0_strided, stride);
  ops_write_plane_group_hdf5({{2, block0np2 / 2}}, slice_name0, {{x0_B0_strided, x1_B0_strided}});
}

void ops_write_plane_group_strided(char *slice_name0, int stride[], ops_dat rho_B0, ops_dat rhou0_B0, ops_dat rhou1_B0,
                                   ops_dat rhou2_B0, ops_dat rhoE_B0, ops_dat WENO_filter_B0) {

  double cpu_start0, elapsed_start0;
  ops_timers(&cpu_start0, &elapsed_start0);

  copy_to_strided_dat(rho_B0, rho_B0_strided, stride);
  copy_to_strided_dat(rhou0_B0, rhou0_B0_strided, stride);
  copy_to_strided_dat(rhou1_B0, rhou1_B0_strided, stride);
  copy_to_strided_dat(rhou2_B0, rhou2_B0_strided, stride);
  copy_to_strided_dat(rhoE_B0, rhoE_B0_strided, stride);
  copy_to_strided_dat(WENO_filter_B0, WENO_filter_B0_strided, stride);

  ops_write_plane_group_hdf5({{2, block0np2 / 2}}, slice_name0,
                             {{rho_B0_strided, rhou0_B0_strided, rhou1_B0_strided, rhou2_B0_strided, rhoE_B0_strided,
                               WENO_filter_B0_strided}});

  double cpu_end0, elapsed_end0;
  ops_timers(&cpu_end0, &elapsed_end0);
  ops_printf("-----------------------------------------\n");
  ops_printf("Time to write strided slice HDF5 file %s: %lf\n", slice_name0, elapsed_end0 - elapsed_start0);
  ops_printf("-----------------------------------------\n");
}
