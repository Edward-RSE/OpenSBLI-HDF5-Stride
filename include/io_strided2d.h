/**
 * @brief Convert between 2D indices and linear indices
 *
 * @param[in] i The i index
 * @param[in] j The j index
 * @param[in] size The size of each dimension
 *
 * @returns The linear index
 *
 */
size_t index2d(size_t i, size_t j, int *size) { return i + j * size[0]; }

/**
 * @brief Creates a strided ops_dat, without halo cells, from an existing
 *        ops_dat.
 *
 * @param[in] dat The ops_dat to create a strided ops_dat from
 * @param[in] stride_i The stride in the i direction
 * @param[in] stride_j The stride in the j direction
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
ops_dat create_strided_ops_dat_without_halo_cells(ops_dat dat, size_t stride_i, size_t stride_j) {
  /*
   * Determine the size of the strided dat halos and the actual size of the
   * dataset. In this case, we do not have any halo cells.
   */
  int strided_d_p[] = {0, 0};
  int strided_d_m[] = {0, 0};
  int strided_size[] = {
      (int)ceil(block0np0 / stride_i),
      (int)ceil(block0np1 / stride_j),
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
  int dat_size[] = {dat->size[0], dat->size[1]};
  int dat_start[] = {(int)abs(dat->d_m[0]), (int)abs(dat->d_m[1])};
  int dat_stop[] = {dat_size[0] - dat->d_p[0], dat_size[1] - dat->d_p[1]};

  /*
   * Loop over each dimension, with each loop variable being incremented by the
   * stride for that dimension.
   */
  for (int j = dat_start[1]; j < dat_stop[1]; j += stride_j) {
    for (int i = dat_start[0]; i < dat_stop[0]; i += stride_i) {
      /*
       * Using `index3d` we can calculate the index in the flat 1d buffers
       * given the 3d coordinates. The coordinate transformation for the
       * strided index is rather complex and ugly to look at. The index is
       * multiplied by the element size of ops_dat, because ops_dat->data is
       * a char*, so we are working with bytes rather than structured data.
       */
      const size_t index_dat = index2d(i, j, dat->size) * dat->elem_size;
      const size_t index_stride =
          index2d((i - dat_start[0]) / stride_i, (j - dat_start[1]) / stride_j, strided_dat->size) *
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

/**
 * @brief Creates a strided ops_dat, including halo cells, from an existing
 *        ops_dat.
 *
 * @param[in] dat The ops_dat to create a strided ops_dat from
 * @param[in] stride_i The stride in the i direction
 * @param[in] stride_j The stride in the j direction
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
ops_dat create_strided_ops_dat_with_halo_cells(ops_dat dat, size_t stride_i, size_t stride_j) {
  /*
   * Determine the size of the strided dat halos and the actual size of the
   * dataset. This is taken from the original ops_dat, but it could be hardcoded
   * during code generation since we already know the sizes of the halos and
   * the blocks.
   */
  int strided_d_p[] = {
      dat->d_p[0] - dat->x_pad,
      dat->d_p[1],
  };
  int strided_d_m[] = {
      dat->d_m[0],
      dat->d_m[1],
  };
  int strided_size[] = {
      (int)ceil(block0np0 / stride_i),
      (int)ceil(block0np1 / stride_j),
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
  int dat_size[] = {dat->size[0], dat->size[1]};
  int dat_start[] = {(int)abs(dat->d_m[0]), (int)abs(dat->d_m[1])};
  int dat_stop[] = {dat_size[0] - dat->d_p[0], dat_size[1] - dat->d_p[1]};

  /*
   * We also need to consider the offset in the strided buffer. This is because
   * we need to worry about the halo cells in the negative direction, d_m,
   * otherwise the final output will be missing data and be in the wrong place.
   */
  int offset_in_strided_buf[] = {(int)abs(strided_dat->d_m[0]), (int)abs(strided_dat->d_m[1])};

  /*
   * Loop over each dimension, with each loop variable being incremented by the
   * stride for that dimension.
   */
  for (int j = dat_start[1]; j < dat_stop[1]; j += stride_j) {
    for (int i = dat_start[0]; i < dat_stop[0]; i += stride_i) {
      /*
       * Using `index3d` we can calculate the index in the flat 1d buffers
       * given the 3d coordinates. The coordinate transformation for the
       * strided index is rather complex and ugly to look at. The index is
       * multiplied by the element size of ops_dat, because ops_dat->data is
       * a char*, so we are working with bytes rather than structured data.
       */
      const size_t index_dat = index2d(i, j, dat->size) * dat->elem_size;
      const size_t index_stride = index2d((i - dat_start[0]) / stride_i + offset_in_strided_buf[0],
                                          (j - dat_start[1]) / stride_j + offset_in_strided_buf[1], strided_dat->size) *
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
