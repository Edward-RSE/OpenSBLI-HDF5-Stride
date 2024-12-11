# Implementation Notes

This document goes over the basic implementation for strided data output in OpenSBLI.

The implementation developed is "standalone," in the sense that no code changes have been made to OPS or to OpenSBLI.
Existing parts of the OPS API are used to create strided datasets (or multi-grids, as it may be known in OPS lingo)
which are passed to the already existing HDF5 output functions. As such, it is conceptually similar to how OpenSBLI uses
OPS already in the sense that the new code is a parallel kernel, in addition to other functions/variables responsible
for initialisation, bookkeeping and writing data to disk.

The heart of the algorithm is a parallel loop (`ops_par_loop`) which uses a "restricted" stencil to copy data from a
dataset into a smaller dataset, using a strided data access pattern. This smaller dataset is the strided dataset which
can be passed to *any* of the HDF5 functions to write data to disk. This also means that the entire process is
parallelised.

To use OPS's HDF5 functions, we need to use OPS's data structures. Therefore we have to copy data *manually* to create a
strided dataset which is passed to OPS. This comes with some overhead and limitations which will be discussed later, but
it is a **strict** requirement due to how the grid is partitioned and spread amongst MPI ranks by OPS. However, there
are some benefits. Namely that we do not need to maintain our own library of code, which effectively re-implements what
is already in OPS, and it also means that everything is parallelised "for free" (because we don't need to parallelise it
ourselves.

**We are taking this approach because the `stride` parameter in the HDF5 hyperslab implementation in OPS results in
blank points, rather than those points not being written.**

## Datasets and Stencils

Each variable which is going to be put output with this implementation needs a `ops_dat` dataset which is the size of
the grid taking into account the stride requested. These datasets are used in conjunction with a "restricted" stencil
which enables strided data access patterns in parallel loops/kernels in OPS.

When defining a dataset, the size of the dataset has to be at least `blockXnpY / strideY` big, where X is the block
number and Y is a direction. The dataset also has to be declared with a `stride` parameter, otherwise something happens
with the stencil and HDF5 output routines which results in garbage output. Unfortunately when declaring a strided
dataset with use with a slice, we can't have the slice dimension be set to 1 because the parallel loop to copy data
(which will be discussed later) *has* to match up with the stride of the stencil and the limits of the parallel loops.
If they don't, the program will either crash or garage will be the output. This does mean that there is some redundant
memory allocated when using the slice or slab output functions in conjunction with this striding implementation.

Below is an example of declaring a strided dataset for a 2D problem. There are no ghost cells defined in the example
below because they are not required. However, ghost cells can be included if you wish.

```c
int strided_size[] = {block0np0 / stride[0], block0np1 / stride[1]};
int strided_base[] = {0, 0};
int strided_d_p[] = {0, 0};
int strided_d_m[] = {0, 0};
double *dummy = NULL;
rho_B0_strided = ops_decl_dat(block, 1, strided_size, strided_base, strided_d_m, strided_d_p, stride, dummy, "double",
                              "rho_B0_strided");
```

Two stencils are required. The first stencil is a "regular" stencil for the point (0, 0). The second stencil is a
restrict stencil, which enables strided data access for parallel kernels. When naming the stencil, it **HAS** to have
the string "RESTRICT" in its name. This is because the OPS translator sees this name and generate different code to
regular stencils.

An example of the two stencils required is shown in the code box below.

```c
int stencil_point[] = {0, 0};
stencil2d_00 = ops_decl_stencil(2, 1, stencil_point, "stencil_00");
stencil2d_restrict_00 = ops_decl_restrict_stencil(2, 1, stencil_point, stride, "stencil_RESTRICT_00");
```

## Parallel copy

Creating the strided dataset requires us to copy data from the original "full fat" dataset, to a smaller dataset which
represents the original dataset but with some stride. The copying is handled by a parallel loop as this is, I believe,
the *only* way to get access to the entire dataset when MPI is used.

```c
int iter_range[] = {0, block0np0 / stride[0], 0, block0np1 / stride[1]};

ops_par_loop(restrict_kernel, "restrict_kernel", block, 2, iter_range,
              ops_arg_dat(original_dat, 1, stencil2d_restrict_00, "double", OPS_READ),
              ops_arg_dat(strided_dat, 1, stencil2d_00, "double", OPS_WRITE),
              ops_arg_idx());
```

The iteration range for the parallel loop is the same size as the strided dataset. The `ops_arg_idx()` is only
technically required for CUDA, but it will compile with warnings if you do not include an index parameter. If you are
hoping to take a "strided slice," then you can limit the iteration range for the slice axis.

```c
int slice_index = block0np2 / 2;
int iter_range[] = {0, block0np0 / stride[0], 0, block0np1, slice_index / stride[1], slice_index / stride[1] + 1};

ops_par_loop(restrict_kernel, "restrict_kernel", block, 2, iter_range,
              ops_arg_dat(original_dat, 1, stencil2d_restrict_00, "double", OPS_READ),
              ops_arg_dat(strided_dat, 1, stencil2d_00, "double", OPS_WRITE),
              ops_arg_idx());
```

In the above example, the parallel loop will only copy data from the slice plane. But, importantly, the dataset still
has to be of `size[] = {block0np0 / stride[0], block0np1 / stride[1], block0np2 / stride[2]}` so the iteration range of
the loop and restricted stencil match up.

The kernel for the parallel loop is very simple, we just set the element of the strided dataset to the corresponding
value in the original dataset. We do not need to worry about indexing to the correct element, e.g. (0 + stride), because
the restrict stencil does that for us.

```c
void restrict_kernel(const ACC<double> &original_dat, ACC<double> &strided_dat, const int *idx) {
  strided_dat(0, 0) = original_dat(0, 0);
}
```

## Single precision data

For further savings, it is possible to use a static cast in a kernel to convert from double to single precision, i.e.
from `double` to `float`. The key changes required to do this are to declare the strided `ops_dat` with a data type
of `float` and to update the `ops_par_loop` and arguments for `restrict_kernel` to use `float` instead of `double`.

In the example below, the elements from `original_dat` are statically cast from `double` to `float` in the kernel, with
`strided_dat` now being of type `ACC<float>`. The corresponding argument in `ops_par_loop` has also been updated to
use `"float"` instead of `"double"` in `ops_arg_dat`.

```c
void restrict_kernel(const ACC<double> &original_dat, ACC<float> &strided_dat, const int *idx) {
  strided_dat(0, 0) = (float)original_dat(0, 0);
}

ops_par_loop(restrict_kernel, "restrict_kernel", block, 2, iter_range,
              ops_arg_dat(original_dat, 1, stencil2d_restrict_00, "double", OPS_READ),
              ops_arg_dat(strided_dat, 1, stencil2d_00, "float", OPS_WRITE),
              ops_arg_idx());
```

## How to use a strided dataset with the HDF5 API

Since the strided dataset is an `ops_dat`, it can be used with any of the OPS HDF5 functions. When using the slab and
slice functions, you need to be aware of the slice/slab range as it has to match the dimensions of the strided dataset
rather than the original.

```c
/* Write the entire dataset */
ops_fetch_dat_hdf5_file(rho_B0_strided, name);
/*  Write a slice of the dataset */
ops_write_plane_group_hdf5({{2, (block0np2 / stride) / 2 }}, name, {{rho_B0_strided}});
/* And so on... */
```

## Usage

The implementation is currently in a file named `io_stride.cpp` with two user-facing functions:

- `HDF5_IO_Init_0_opensbliblock00_strided`: initialises the datasets and stencils for strided output
- `HDF5_IO_Write_0_opensbliblock00_strided`: populates the datasets and writes them to disk

The other functions and variables in `io_stride.cpp` are designed to be hidden away from the user. However, there is no
reason why they should be other than to hide away complication in the main function. The init function has to be called
before `ops_partition`, otherwise OPS will not be able to distribute data across MPI ranks and will exit.

```c
/* Initialisation only requires the block, block dimensions and the stride parameters.
   The strided datasets are in `io_stride.cpp` and don't need to be declared in scope of
   the main function in this implementation */

int output_stride = {2, 2};
HDF5_IO_Init_0_opensbliblock00_strided(opensbliblock00, block0np0, block0np1, output_stride);

/* Other initialisation code goes here */

ops_partition("");

/* To write strided output, just call this function with the correct arguments.
   rho_B0 is the ops_dat for the density field which will be transformed into a
   strided dataset and written to disk */

HDF5_IO_Write_0_opensbliblock00_strided(opensbliblock00, block0np0, block0np1, output_stride, &rho_B0)
```

## Limitations

1. The `ops_dat` for the strided dataset has to be large enough for the iteration range of the parallel loop w.r.t. the
   size of the original dataset. Although you may be able to work around this with a stride like `{2, 2, block0np2 /
   2}`. For example, for a strided slice with the slice happening on the z axis (direction 2), you can't have an ops_dat
   with size `{block0np0 / stride, block0np1 / stride, 1}`. The parallel loop or the HDF5 output will not be happy,
   unless you juggle around a weird stencil and have separate ranges for the output and parallel loop.
2. In MPI modes, have to declare all the `ops_dat`'s before `ops_partition()`. So you cannot dynamically allocate and
   free the memory for the strided `ops_dat`'s. This increases memory requirements. Of course, this is not so bad if the
   stride is large enough.
