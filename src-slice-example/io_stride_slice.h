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

/**
 * @brief Creates a strided ops_dat, without halo cells, from an existing
 *        ops_dat.
 *
 * @param[in] dat The ops_dat to create a strided ops_dat from
 * @param[in] stride_i The stride in the i direction
 * @param[in] stride_j The stride in the j direction
 * @param[in] stride_k The stride in the k direction
 *
 * @returns The strided ops_dat
 *
 * @details
 *
 * This function creates a strided ops_dat by taking an existing ops_dat and
 * creating a new ops_dat with the specified stride in each dimension. There
 * are no halo cells, because it simplifies copying data between buffers.
 *
 */
ops_dat create_strided_ops_dat_without_halo_cells(ops_dat dat, int stride[3]) {
  /*
   * Determine the size of the strided dat halos and the actual size of the
   * dataset. In this case, we do not have any halo cells.
   */
  int strided_d_p[] = {0, 0, 0};
  int strided_d_m[] = {0, 0, 0};
  int strided_size[] = {
      (int)ceil((dat->size[0] - dat->d_p[0] + dat->d_m[0]) / stride[0]),
      (int)ceil((dat->size[1] - dat->d_p[1] + dat->d_m[1]) / stride[1]),
      (int)ceil(dat->size[2] - dat->d_p[2] + dat->d_m[2]),
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
  double *dummy = NULL;
  ops_dat strided_dat = ops_decl_dat(dat->block, dat->dim, strided_size, dat->base, strided_d_m, strided_d_p, dummy,
                                     dat->type, dat->name);

  /*
   * Loop variables to control the start and stop of the loop over each
   * dimension. These values are based on the size of the original dataset. The
   * loops start at the first non-halo cell and stop at the last non-halo cell.
   * The loop ranges are over the entire dataset, dat_size.
   */
  int dat_size[] = {dat->size[0], dat->size[1], dat->size[2]};
  int dat_start[] = {(int)abs(dat->d_m[0]), (int)abs(dat->d_m[1]), (int)abs(dat->d_m[2])};
  int dat_stop[] = {dat_size[0] - dat->d_p[0], dat_size[1] - dat->d_p[1], dat_size[2] - dat->d_p[2]};

  /*
   * Loop over each dimension, with each loop variable being incremented by the
   * stride for that dimension.
   */

  const int k = (dat->size[2] - dat->d_p[2] + dat->d_m[2]) / 2;

  for (int j = dat_start[1]; j < dat_stop[1]; j += stride[1]) {
    for (int i = dat_start[0]; i < dat_stop[0]; i += stride[0]) {
      /*
       * Using `index3d` we can calculate the index in the flat 1d buffers
       * given the 3d coordinates. The coordinate transformation for the
       * strided index is rather complex and ugly to look at. The index is
       * multiplied by the element size of ops_dat, because ops_dat->data is
       * a char*, so we are working with bytes rather than structured data.
       */
      size_t index_dat = index3d(i, j, k, dat->size) * dat->elem_size;
      size_t index_stride =
          index3d((i - dat_start[0]) / stride[0], (j - dat_start[1]) / stride[1], k, strided_dat->size) *
          strided_dat->elem_size;

      /*
       * Copy the data from the original dataset into the strided buffer. This
       * has to be done using memcpy because ops_dat->data is a char* so we
       * are working in bytes, rather than structured data.
       */
      memcpy(&strided_dat->data[index_stride], &dat->data[index_dat], strided_dat->elem_size);
    }
  }

  return strided_dat;
}

void ops_write_plane_group_strided_coords(char *slice_name0, ops_dat x0_B0, ops_dat x1_B0) {
  int stride[] = {2, 2, 1};
  ops_dat x0_B0_strided = create_strided_ops_dat_without_halo_cells(x0_B0, stride);
  ops_dat x1_B0_strided = create_strided_ops_dat_without_halo_cells(x1_B0, stride);
  ops_write_plane_group_hdf5({{2, block0np2 / 2}}, slice_name0, {{x0_B0_strided, x1_B0_strided}});
}

void ops_write_plane_group_strided(char *slice_name0, ops_dat rho_B0, ops_dat rhou0_B0, ops_dat rhou1_B0,
                                   ops_dat rhou2_B0, ops_dat rhoE_B0, ops_dat WENO_filter_B0) {

  int stride[] = {2, 2, 1};

  ops_dat rho_B0_strided = create_strided_ops_dat_without_halo_cells(rho_B0, stride);
  ops_dat rhou0_B0_strided = create_strided_ops_dat_without_halo_cells(rhou0_B0, stride);
  ops_dat rhou1_B0_strided = create_strided_ops_dat_without_halo_cells(rhou1_B0, stride);
  ops_dat rhou2_B0_strided = create_strided_ops_dat_without_halo_cells(rhou2_B0, stride);
  ops_dat rhoE_B0_strided = create_strided_ops_dat_without_halo_cells(rhoE_B0, stride);
  ops_dat WENO_filter_B0_strided = create_strided_ops_dat_without_halo_cells(WENO_filter_B0, stride);

  ops_write_plane_group_hdf5({{2, block0np2 / 2}}, slice_name0,
                             {{rho_B0_strided, rhou0_B0_strided, rhou1_B0_strided, rhou2_B0_strided, rhoE_B0_strided,
                               WENO_filter_B0_strided}});
}
