"""Compare reference data against output data for a wave equation OpenSBLI
simulation.
"""

import re
import argparse
import sys
from pathlib import Path

import h5py
import numpy
from matplotlib import pyplot


def main(args):
    try:
        original_dat = h5py.File(args.original_data, "r")
    except IOError:
        print(f"Unable to open: {args.original_data}")
        sys.exit(1)
    try:
        strided_dat = h5py.File(args.new_data, "r")
    except IOError:
        print(f"Unable to open: {args.new_data}")
        sys.exit(1)

    variable = "rho_B0"
    n_iter = int(re.search(r"\d+", args.original_data).group())
    halo_n = 5
    halo_p = 5

    variable_original = numpy.array(
        original_dat["opensbliblock00"][variable][
            halo_n:-halo_p, halo_n:-halo_p, halo_n:-halo_p
        ]
    )
    variable_slab_2d = numpy.array(
        strided_dat["opensbliblock00"][f"{n_iter}"][f"{variable}"][:, :, :]
    )

    offset0 = 40
    offset1 = 40
    offset2 = 10
    block0np0 = 180
    block0np1 = 175
    block0np2 = 150

    start0 = int((block0np2 / 2) - offset2)
    end0 = int((block0np2 / 2) + offset2)
    start1 = int((block0np1 / 2) - offset1)
    end1 = int((block0np1 / 2) + offset1)
    start2 = int((block0np0 / 2) - offset0)
    end2 = int((block0np0 / 2) + offset0)

    # start0 = 0
    # end0 = block0np2
    # start1 = 0
    # end1 = block0np1
    # start2 = 0
    # end2 = block0np0

    variable_original_slice = variable_original[start0:end0, start1:end1, start2:end2]
    variable_original_slice = variable_original_slice[
        :: args.stride_i, :: args.stride_j, :: args.stride_k
    ]

    print("Slab range:", start0, end0, start1, end1, start2, end2)
    print(
        args.original_data,
        "before slice",
        variable_original.shape,
        "after slice",
        variable_original_slice.shape,
    )
    print(args.new_data, variable_slab_2d.shape)

    for k in range(variable_slab_2d.shape[0]):
        fig, ax = pyplot.subplots(1, 5, figsize=(12, 3))
        # Original grid, no modifications
        ax[0].imshow(variable_original[k, :, :])
        ax[0].set_title("Original grid")
        # Original grid, but used slicing to get the slab
        ax[1].imshow(variable_original_slice[k, :, :])
        ax[1].set_title("Original slab [sliced in numpy]")
        # The output from the slab HDF5
        ax[2].imshow(variable_slab_2d[k, :, :])
        ax[2].set_title(
            f"HDF5 slab [stride = \n{args.stride_i}, {args.stride_j}, {args.stride_k}]"
        )
        # Absolute difference between the original grid, sliced in numpy, and
        # the output form the slab HDF5
        diff = variable_original_slice[k, :, :] - variable_slab_2d[k, :, :]
        if numpy.any(diff > 0):
            print(f"Difference for k = {k}")
        ax[3].imshow(diff)
        ax[3].set_title("Difference = slab -\n strided slab")
        # Absolute difference between the original grid and the grid sliced in
        # numpy to create a slab
        diff = (
            variable_original[k + start0, start1:end1, start2:end2]
            - variable_original_slice[k, :, :]
        )
        ax[4].imshow(diff)
        ax[4].set_title("Difference = original -\n original slab")
        # Rest of the figure
        fig.suptitle(f"Slice index: {k}")
        fig.tight_layout()
        path = Path(__file__).parent
        fig.savefig(f"{path}/fig/slab-example_k={k}.png")

        pyplot.close()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("original_data")
    ap.add_argument("new_data")
    ap.add_argument("stride_i", default=1, type=int)
    ap.add_argument("stride_j", default=1, type=int)
    ap.add_argument("stride_k", default=1, type=int)
    args = ap.parse_args()
    main(args)
