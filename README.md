# Data Compression via Strided Output in OpenSBLI

HPC RSE project at the University of Southampton. The aim is to develop a method
for reducing the size of output HDF5 files from an OpenSBLI simulation by
using striding to writing to disk at specified intervals of the grid -- such as
every second, fourth or nth grid point -- instead of saving data from every
grid point.

The implementation is described in `docs/implementation.md` and the code in each
directory is in `io_stride.cpp`.

## Building

The code in this repository does not require OpenSBLI, as we will be prototyping
for pre-generated code. To run (and develop) these example, you will need
installed on your system:

- OPS
- Parallel HDF5
- A parallelisation framework, e.g. MPI or CUDA

Each example should have a `build.sh` build script which can be used to
build the examples, note that you will be required to update the scripts to
point to your OPS and HDF5 installation locations.
