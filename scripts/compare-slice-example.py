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
    block0np2 = 150
    halo_n = 5
    halo_p = 5

    original_dat_path = Path(args.original_data)
    original_dat_iter = int(original_dat_path.stem.split("_")[-1])
    variable_original = numpy.array(
        original_dat["opensbliblock00"][variable][
            halo_n:-halo_p, halo_n:-halo_p, halo_n:-halo_p
        ]
    )
    variable_strided = numpy.array(
        strided_dat["opensbliblock00"][f"{original_dat_iter}"][variable][:, :]
    )

    slice_index = int((block0np2) / 2)
    variable_original_slice = variable_original[slice_index, :, :]
    variable_original_slice = variable_original_slice[
        :: args.stride_j, :: args.stride_i
    ]

    diff = variable_original_slice - variable_strided

    fig, ax = pyplot.subplots(1, 3, figsize=(12, 5))
    ax[0].imshow(variable_original[slice_index, :, :])
    ax[1].imshow(variable_strided)
    im = ax[2].imshow(diff)
    fig.colorbar(im, ax=ax[2])
    ax[0].set_title("Original slice")
    ax[1].set_title("Strided [slice with striding]")
    ax[2].set_title("Absolute difference")
    fig.suptitle(f"Slice index: {slice_index}")
    fig.tight_layout()

    pyplot.show()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("original_data")
    ap.add_argument("new_data")
    ap.add_argument("stride_i", default=1, type=int)
    ap.add_argument("stride_j", default=1, type=int)
    ap.add_argument("stride_k", default=1, type=int)
    args = ap.parse_args()
    main(args)
