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
ops_dat create_strided_ops_dat_without_halo_cells(ops_dat dat, size_t stride_i, size_t stride_j, size_t stride_k) {
  /*
   * Determine the size of the strided dat halos and the actual size of the
   * dataset. In this case, we do not have any halo cells.
   */
  int strided_d_p[] = {0, 0, 0};
  int strided_d_m[] = {0, 0, 0};
  int strided_size[] = {
      (int)ceil(block0np0 / stride_i),
      (int)ceil(block0np1 / stride_j),
      (int)ceil(block0np2 / stride_k),
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
  for (int k = dat_start[2]; k < dat_stop[2]; k += stride_k) {
    for (int j = dat_start[1]; j < dat_stop[1]; j += stride_j) {
      for (int i = dat_start[0]; i < dat_stop[0]; i += stride_i) {
        /*
         * Using `index3d` we can calculate the index in the flat 1d buffers
         * given the 3d coordinates. The coordinate transformation for the
         * strided index is rather complex and ugly to look at. The index is
         * multiplied by the element size of ops_dat, because ops_dat->data is
         * a char*, so we are working with bytes rather than structured data.
         */
        const size_t index_dat = index3d(i, j, k, dat->size) * dat->elem_size;
        const size_t index_stride = index3d((i - dat_start[0]) / stride_i, (j - dat_start[1]) / stride_j,
                                            (k - dat_start[2]) / stride_k, strided_dat->size) *
                                    strided_dat->elem_size;
        /*
         * Copy the data from the original dataset into the strided buffer. This
         * has to be done using memcpy because ops_dat->data is a char* so we
         * are working in bytes, rather than structured data.
         */
        memcpy(&strided_dat->data[index_stride], &dat->data[index_dat], strided_dat->elem_size);
      }
    }
  }

  return strided_dat;
}

/**
 * @brief Creates a strided ops_dat, including halo cells, from an existing
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
 * creating a new ops_dat with the specified stride in each dimension. The same
 * number of halo cells are used in each dimension, taken from the original
 * ops_dat, but will contain uninitialized (or zero'd) data. This is done for
 * consistency and to simplify plotting, as you do not need to worry about
 * different starting indices and padding when processing the data.
 *
 */
ops_dat create_strided_ops_dat_with_halo_cells(ops_dat dat, size_t stride_i, size_t stride_j, size_t stride_k) {
  /*
   * Determine the size of the strided dat halos and the actual size of the
   * dataset. This is taken from the original ops_dat, but it could be hardcoded
   * during code generation since we already know the sizes of the halos and
   * the blocks.
   */
  int strided_d_p[] = {
      dat->d_p[0] - dat->x_pad,
      dat->d_p[1],
      dat->d_p[2],
  };
  int strided_d_m[] = {
      dat->d_m[0],
      dat->d_m[1],
      dat->d_m[2],
  };
  int strided_size[] = {
      (int)ceil(block0np0 / stride_i),
      (int)ceil(block0np1 / stride_j),
      (int)ceil(block0np2 / stride_k),
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
   * We also need to consider the offset in the strided buffer. This is because
   * we need to worry about the halo cells in the negative direction, d_m,
   * otherwise the final output will be missing data and be in the wrong place.
   */
  int offset_in_strided_buf[] = {(int)abs(strided_dat->d_m[0]), (int)abs(strided_dat->d_m[1]),
                                 (int)abs(strided_dat->d_m[2])};

  /*
   * Loop over each dimension, with each loop variable being incremented by the
   * stride for that dimension.
   */
  for (int k = dat_start[2]; k < dat_stop[2]; k += stride_k) {
    for (int j = dat_start[1]; j < dat_stop[1]; j += stride_j) {
      for (int i = dat_start[0]; i < dat_stop[0]; i += stride_i) {
        /*
         * Using `index3d` we can calculate the index in the flat 1d buffers
         * given the 3d coordinates. The coordinate transformation for the
         * strided index is rather complex and ugly to look at. The index is
         * multiplied by the element size of ops_dat, because ops_dat->data is
         * a char*, so we are working with bytes rather than structured data.
         */
        const size_t index_dat = index3d(i, j, k, dat->size) * dat->elem_size;
        const size_t index_stride =
            index3d((i - dat_start[0]) / stride_i + offset_in_strided_buf[0],
                    (j - dat_start[1]) / stride_j + offset_in_strided_buf[1],
                    (k - dat_start[2]) / stride_k + offset_in_strided_buf[2], strided_dat->size) *
            strided_dat->elem_size;
        /*
         * Copy the data from the original dataset into the strided buffer. This
         * has to be done using memcpy because ops_dat->data is a char* so we
         * are working in bytes, rather than structured data.
         */
        memcpy(&strided_dat->data[index_stride], &dat->data[index_dat], strided_dat->elem_size);
      }
    }
  }

  return strided_dat;
}
