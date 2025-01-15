# Performance

The purpose of adding a stride parameter to data output is to reduce file sizes.

## Methodology

Grid size is 650 x 610 x 520, which is ~206 million grid points. Just 15 iterations because we don't care about how
long the simulation takes, only interested in the time spent in I/O and the size of the final output. We use a stride of
{5, 5}, which means a slice in the z plane will have dimensions 130 x 101.

Only writing out one variable, the density field `rho_B0`, in multiple ways: 1) writing the entire grid to disk, 2)
writing a slice to disk, and 3) writing a strided slice to disk. We are using a parallel filesystem. Always writing to
disk in double precision, not going to worry about single.

Using two processes and two A100 GPUs, so the I/O is parallel in this case.

## Benchmark

| Method                   | Time (ms) | File size (MB) |
| ------------------------ | --------- | -------------- |
| Regular                  | 3455      | 1700           |
| Slice                    | 92        | 6              |
| Slice with stride {5, 5} | 178       | 0.25           |

The times are averaged over 15 events. The file size of regular is skewed because it has to include the block, which
must be adding in lots of additional data.

However, I have also found setups where the strided slice output is faster than the regular slice output. For example,
it is roughly 3-10x faster in a setup where the output frequency is decreased (e.g. write every 200 iterations rather
than every iteration):

| Method                   | Time (ms) |
| ------------------------ | --------- |
| Slice                    | 166       |
| Slice with stride {5, 5} | 40        |

The simulation ran on the same hardware and setup, and the results are an average over 35 calls. I'm unsure why the
regular slice took longer this time.

## Optimisation

The underlying issue for the strided output being slower is because it copies data from one `ops_dat` to
another which the regular output methods don't do. So we are doing extra work. There's no way around this, unless
changes are made to the OPS HDF5 API. It would require either changing `fetch_loop_slab()` to include a `stride`
argument and accounting for it in the various loops, or by unifying the HDF5 to use HDF5 Hyperslabs which include a
stride parameter. Hyperslabs already are used in OPS, but only seem to be present in the MPI API.

## Old benchmark example

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