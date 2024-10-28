"""Compare reference data against output data for a wave equation OpenSBLI
simulation.
"""

import argparse
import sys
from pathlib import Path

import h5py
import numpy
from IPython import embed
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

    halo_n = 5
    halo_p = 5
    rho_original = numpy.array(
        original_dat["opensbliblock00"]["rho_B0"][
            halo_n:-halo_p, halo_n:-halo_p, halo_n:-halo_p
        ]
    )

    original_dat_path = Path(args.original_data)
    dat_iter = int(original_dat_path.stem.split("_")[-1])
    rho_strided = numpy.array(
        strided_dat["opensbliblock00"][f"{dat_iter}"]["rho_B0"][:, :]
    )

    block0np2 = 64
    # slice_index = int((block0np2 / args.stride_k) / 2)
    slice_index = int(block0np2 / 2)
    rho_original_slice = rho_original[:, :, slice_index]
    rho_original_slice_strided = rho_original_slice[:: args.stride_i, :: args.stride_j]

    print(slice_index)
    print(rho_strided.shape, rho_original_slice_strided.shape)

    diff = rho_original_slice_strided - rho_strided

    fig, ax = pyplot.subplots(1, 3, figsize=(12, 5))

    im0 = ax[0].imshow(rho_original_slice_strided)
    im1 = ax[1].imshow(rho_strided)
    im2 = ax[2].imshow(diff)

    # fig.colorbar(im0, ax=ax[0])
    # fig.colorbar(im1, ax=ax[1])
    # fig.colorbar(im2, ax=ax[2])

    ax[0].set_title("Original [manually sliced and strided]")
    ax[1].set_title("Strided [slice with striding]")
    ax[2].set_title("Absolute difference")

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
