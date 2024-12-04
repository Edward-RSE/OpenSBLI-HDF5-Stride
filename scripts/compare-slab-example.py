"""Compare reference data against output data for a wave equation OpenSBLI
simulation.
"""

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
    n_iter = 4
    halo_n = 5
    halo_p = 5

    variable_original = numpy.array(
        original_dat["opensbliblock00"][variable][
            # :, :, :
            halo_n:-halo_p, halo_n:-halo_p, halo_n:-halo_p
        ]
    )
    variable_slab_2d = numpy.array(
        strided_dat["opensbliblock00"][f"{n_iter}"][f"{variable}"][
            :, :, :
            # halo_n:-halo_p, halo_n:-halo_p, halo_n:-halo_p
        ]
    )

    block0np0 = 180
    block0np1 = 175
    block0np2 = 150
    variable_original_slice = variable_original[0:block0np0, 0:block0np1, 0:block0np2]
    variable_original_slice = variable_original

    print(args.original_data, variable_original_slice.shape)
    print(args.new_data, variable_slab_2d.shape)

    for k in range(variable_slab_2d.shape[2]):
        fig, ax = pyplot.subplots(1, 3, figsize=(12, 5))
        ax[0].imshow(variable_original[:, :, k])
        ax[1].imshow(variable_slab_2d[:, :, k])
        diff = variable_original_slice[:, :, k] - variable_slab_2d[:, :, k]
        if numpy.any(diff > 0):
            print(f"Difference for k = {k}")
        im = ax[2].imshow(diff)
        fig.colorbar(im, ax=ax[2])
        ax[0].set_title("Original slice")
        ax[1].set_title("Strided [slice with striding]")
        ax[2].set_title("Absolute difference")
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
