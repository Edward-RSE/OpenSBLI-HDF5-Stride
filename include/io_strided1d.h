/*
 * Takes an input dataset to create a new strided dataset
 * In this example, we create a new ops_dat without any halo cells, which
 * significantly simplifies copying data between buffers.
 */
ops_dat create_strided_ops_dat(ops_dat dat, size_t stride_i) {

  int d_p[] = {0};
  int d_m[] = {0};
  int new_size[] = {ceil(block0np0 / stride_i)};

  double *_dummy = NULL; /* Note that this has to be the type you want the ops_dat to be */
  ops_dat new_dat = ops_decl_dat(dat->block, dat->dim, new_size, dat->base, d_m, d_p, _dummy, dat->type, dat->name);

  int size[] = {dat->size[0]};
  int offset[] = {abs(dat->d_m[0])};

  for (int i = offset[0]; i < size[0] - dat->d_p[0]; i += stride_i) {
    const size_t index_dat = i * dat->elem_size;
    const size_t index_new = (i - offset[0]) * new_dat->elem_size;
    memcpy(new_dat->data + index_new, dat->data + index_dat, new_dat->elem_size);
  }

  return new_dat;
}
