/**
 * @brief Convert between 1D indices and linear indices
 *
 * @param[in] i The i index
 * @param[in] size The size of each dimension (unused)
 *
 * @returns The linear index
 *
 */
size_t index1d(size_t i, int *size) { return i; }

/**
 * @brief Creates a strided ops_dat, without halo cells, from an existing
 *        ops_dat.
 *
 * @param[in] dat The ops_dat to create a strided ops_dat from
 * @param[in] stride_i The stride in the i direction
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
ops_dat create_strided_ops_dat_without_halo_cells(ops_dat dat, size_t stride_i) {
  /*
   * Determine the size of the strided dat halos and the actual size of the
   * dataset. In this case, we do not have any halo cells.
   */
  int strided_d_p[] = {0};
  int strided_d_m[] = {0};
  int strided_size[] = {(int)ceil(block0np0 / stride_i)};

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
  int dat_size[] = {dat->size[0]};
  int dat_start[] = {(int)abs(dat->d_m[0])};
  int dat_stop[] = {dat_size[0] - dat->d_p[0]};

  /*
   * Loop over each dimension, with each loop variable being incremented by the
   * stride for that dimension.
   */
  for (int i = dat_start[0]; i < dat_stop[0]; i += stride_i) {
    /*
     * Using `index3d` we can calculate the index in the flat 1d buffers
     * given the 3d coordinates. The coordinate transformation for the
     * strided index is rather complex and ugly to look at. The index is
     * multiplied by the element size of ops_dat, because ops_dat->data is
     * a char*, so we are working with bytes rather than structured data.
     */
    const size_t index_dat = index1d(i, dat->size) * dat->elem_size;
    const size_t index_stride = index1d((i - dat_start[0]) / stride_i, strided_dat->size) * strided_dat->elem_size;
    /*
     * Copy the data from the original dataset into the strided buffer. This
     * has to be done using memcpy because ops_dat->data is a char* so we
     * are working in bytes, rather than structured data.
     */
    memcpy(&strided_dat->data[index_stride], &dat->data[index_dat], strided_dat->elem_size);
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
ops_dat create_strided_ops_dat_with_halo_cells(ops_dat dat, size_t stride_i) {
  /*
   * Determine the size of the strided dat halos and the actual size of the
   * dataset. This is taken from the original ops_dat, but it could be hardcoded
   * during code generation since we already know the sizes of the halos and
   * the blocks.
   */
  int strided_d_p[] = {dat->d_p[0] - dat->x_pad};
  int strided_d_m[] = {dat->d_m[0]};
  int strided_size[] = {(int)ceil(block0np0 / stride_i)};

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
  int dat_size[] = {dat->size[0]};
  int dat_start[] = {(int)abs(dat->d_m[0])};
  int dat_stop[] = {dat_size[0] - dat->d_p[0]};

  /*
   * We also need to consider the offset in the strided buffer. This is because
   * we need to worry about the halo cells in the negative direction, d_m,
   * otherwise the final output will be missing data and be in the wrong place.
   */
  int offset_in_strided_buf[] = {(int)abs(strided_dat->d_m[0])};

  /*
   * Loop over each dimension, with each loop variable being incremented by the
   * stride for that dimension.
   */
  for (int i = dat_start[0]; i < dat_stop[0]; i += stride_i) {
    /*
     * Using `index3d` we can calculate the index in the flat 1d buffers
     * given the 3d coordinates. The coordinate transformation for the
     * strided index is rather complex and ugly to look at. The index is
     * multiplied by the element size of ops_dat, because ops_dat->data is
     * a char*, so we are working with bytes rather than structured data.
     */
    const size_t index_dat = index1d(i, dat->size) * dat->elem_size;
    const size_t index_stride =
        index1d((i - dat_start[0]) / stride_i + offset_in_strided_buf[0], strided_dat->size) * strided_dat->elem_size;
    /*
     * Copy the data from the original dataset into the strided buffer. This
     * has to be done using memcpy because ops_dat->data is a char* so we
     * are working in bytes, rather than structured data.
     */
    memcpy(&strided_dat->data[index_stride], &dat->data[index_dat], strided_dat->elem_size);
  }

  return strided_dat;
}

/**
 * @brief Creates a strided ops_dat, including halo cells, from an existing
 *        ops_dat converting double to single precision.
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
// ops_dat create_strided_ops_dat_with_halo_cells_single_precision(ops_dat dat, size_t stride_i) {
//   /*
//    * Determine the size of the strided dat halos and the actual size of the
//    * dataset. This is taken from the original ops_dat, but it could be hardcoded
//    * during code generation since we already know the sizes of the halos and
//    * the blocks.
//    */
//   int strided_d_p[] = {dat->d_p[0] - dat->x_pad};
//   int strided_d_m[] = {dat->d_m[0]};
//   int strided_size[] = {(int)ceil(block0np0 / stride_i)};

//   /*
//    * Create a new dataset using the sizes above. We need to create a new ops_dat
//    * rather than just a buffer since all the HDF5 functions (which we want to
//    * re-use) expect an ops_dat and not a buffer containing data. Note that we
//    * have a dummy variable. This variable *HAS* to be the type of data which
//    * will be stored in the ops_dat. If it is not, this will allocate an
//    * incorrect amount of data and will also initialise some of the fields, such
//    * as ops_dat->elem_size, incorrectly. Anything wrong there will result in
//    * copy errors and garbage output. We are specifying FLOAT here because we
//    * want to convert from double to float, e.g. from double to single precision.
//    */
//   float *dummy = NULL;
//   ops_dat strided_dat =
//       ops_decl_dat(dat->block, dat->dim, strided_size, dat->base, strided_d_m, strided_d_p, dummy, "float",
//       dat->name);

//   /*
//    * Loop variables to control the start and stop of the loop over each
//    * dimension. These values are based on the size of the original dataset. The
//    * loops start at the first non-halo cell and stop at the last non-halo cell.
//    * The loop ranges are over the entire dataset, dat_size.
//    */
//   int dat_size[] = {dat->size[0]};
//   int dat_start[] = {(int)abs(dat->d_m[0])};
//   int dat_stop[] = {dat_size[0] - dat->d_p[0]};

//   /*
//    * We also need to consider the offset in the strided buffer. This is because
//    * we need to worry about the halo cells in the negative direction, d_m,
//    * otherwise the final output will be missing data and be in the wrong place.
//    */
//   int offset_in_strided_buf[] = {(int)abs(strided_dat->d_m[0])};

//   /*
//    * Loop over each dimension, with each loop variable being incremented by the
//    * stride for that dimension.
//    */
//   for (int i = dat_start[0]; i < dat_stop[0]; i += stride_i) {
//     /*
//      * Using `index3d` we can calculate the index in the flat 1d buffers
//      * given the 3d coordinates. The coordinate transformation for the
//      * strided index is rather complex and ugly to look at. The index is
//      * multiplied by the element size of ops_dat, because ops_dat->data is
//      * a char*, so we are working with bytes rather than structured data.
//      */
//     const size_t index_dat = index1d(i, dat->size) * dat->elem_size;
//     const size_t index_stride =
//         index1d((i - dat_start[0]) / stride_i + offset_in_strided_buf[0], strided_dat->size) *
//         strided_dat->elem_size;

//     double tmp_d;
//     memcpy(&tmp_d, &dat->data[index_dat], dat->elem_size);
//     const float tmp_f = (float)tmp_d;

//     /*
//      * Copy the data from the original dataset into the strided buffer. This
//      * has to be done using memcpy because ops_dat->data is a char* so we
//      * are working in bytes, rather than structured data.
//      */
//     memcpy(&strided_dat->data[index_stride], &tmp_f, strided_dat->elem_size);
//   }

//   return strided_dat;
// }

#include "io_strided_util.h"

ops_dat create_strided_ops_dat_with_halo_cells_single_precision(ops_dat dat, size_t stride_i) {

  int strided_d_p[] = {dat->d_p[0] - dat->x_pad};
  int strided_d_m[] = {dat->d_m[0]};
  int strided_size[] = {(int)ceil(block0np0 / stride_i)};

  ops_dat strided_dat = create_empty_ops_dat<float>(dat->block, dat->dim, strided_size, dat->base, strided_d_m,
                                                    strided_d_p, "float", dat->name);

  int dat_size[] = {dat->size[0]};
  int dat_start[] = {(int)abs(dat->d_m[0])};
  int dat_stop[] = {dat_size[0] - dat->d_p[0]};

  int offset_in_strided_buf[] = {(int)abs(strided_dat->d_m[0])};

  for (int i = dat_start[0]; i < dat_stop[0]; i += stride_i) {

    const size_t index_dat = index1d(i, dat->size) * dat->elem_size;
    const size_t index_stride =
        index1d((i - dat_start[0]) / stride_i + offset_in_strided_buf[0], strided_dat->size) * strided_dat->elem_size;
    const float float_value = convert_precision<double, float>(&dat->data[index_dat], dat->elem_size);
    memcpy(&strided_dat->data[index_stride], &float_value, strided_dat->elem_size);
  }

  return strided_dat;
}
