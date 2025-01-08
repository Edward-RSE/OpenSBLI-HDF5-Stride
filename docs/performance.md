# Performance

The purpose of adding a stride parameter to data output is to reduce file sizes and time spent writing to disk.

## Stride example

## Slice with stride example

This is a test of using the HDF5 slice option with a stride of {2, 2, 1}, where the 1 is the slice plane. The size of
the grid is 180 x 190 x 150, with the slice plane being nz = 75.

| Method            | Time (ms) | File size (MB) |
| ----------------- | --------- | -------------- |
| Regular           | 120       | 42.0           |
| Slice             | 0.529     | 0.484          |
| Slice with stride | 0.485     | 0.125          |

The time for "slice" and "slice with stride" varies from run-to-run. Typically the slide with stride eeked out victory,
but if did sometimes take longer. In general, there is little performance difference between them with a stride of
{2, 2, 1}. The difference is in the file size, which is 4x smaller for "slice with stride."

Let's now try with a grid of 360 x 355 x 300 and  stride of {5, 5, 1}.

| Method            | Time (ms) | File size (MB) |
| ----------------- | --------- | -------------- |
| Regular           | 321       | 320            |
| Slice             | 1.109     | 2              |
| Slice with stride | 2.239     | 0.086          |

Clearly, the data saving is very good. However, it takes a long time to copy the data to the smaller datasets, which is
not ideal. Let's see how this scales by doing the copy in parallel, using a GPU on Iridis X.

| Method            | Time (ms) | File size (MB) |
| ----------------- | --------- | -------------- |
| Regular           | 211       | 320            |
| Slice             | 4.722     | 2              |
| Slice with stride | 8.119     | 0.086          |

Still takes longer. Although, I'm not sure why it takes longer to write a slice to disk on Iridis X, when it takes less
time to write the "regular" output to disk. I presume this is something to do with the parallel file system.

## Optimisation

The underlying issue for the strided output being slower, or comparable, is because it copies data from one `ops_dat` to
another which the regular output methods don't do. So we are doing extra work. There's no way around this, unless
changes are made to the OPS HDF5 API. It would require either changing `fetch_loop_slab()` to include a `stride`
argument and accounting for it in the various loops, or by unifying the HDF5 to use HDF5 Hyperslabs which include a
stride parameter. Hyperslabs already are used in OPS, but only seem to be present in the MPI API.

Note the slice here is using the stride implementation with a stride of {1 , 1, 1} so it would be faster using the
intended way to create a slice.

Regular 36.67 s
slice 0.351 s
slice with stride