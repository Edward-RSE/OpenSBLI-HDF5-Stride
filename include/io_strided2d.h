/*
 * Helper function to convert indices to a flattened index
 */
size_t index2d(size_t i, size_t j, int *size) { return i + j * size[0]; }

/*
 * Takes an input dataset to create a new strided dataset
 * In this example, we create a new ops_dat without any halo cells, which
 * significantly simplifies copying data between buffers.
 */
ops_dat create_strided_ops_dat(ops_dat dat, size_t stride_i, size_t stride_j, ) {
  int d_p[] = {0, 0};
  int d_m[] = {0, 0};
  int new_size[] = {(int)ceil(block0np0 / stride_i), (int)ceil(block0np1 / stride_j)};

  double *_dummy = NULL; /* Note that this has to be the type you want the ops_dat to be */
  ops_dat new_dat = ops_decl_dat(dat->block, dat->dim, new_size, dat->base, d_m, d_p, _dummy, dat->type, dat->name);

  int size[] = {dat->size[0], dat->size[1]};
  int offset[] = {(int)abs(dat->d_m[0]), (int)abs(dat->d_m[1])};

  for (int j = offset[1]; j < size[1] - dat->d_p[1]; j += stride_j) {
    for (int i = offset[0]; i < size[0] - dat->d_p[0]; i += stride_i) {
      const size_t index_dat = index2d(i, j, dat->size) * dat->elem_size;
      const size_t index_new =
          index3d((i - offset[0]) / stride_i, (j - offset[1]) / stride_j, new_dat->size) * new_dat->elem_size;
      memcpy(new_dat->data + index_new, dat->data + index_dat, new_dat->elem_size);
    }
  }

  return new_dat;
}
