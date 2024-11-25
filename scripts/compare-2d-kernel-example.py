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

    halo_n = 5
    halo_p = 5
    rho_original = numpy.array(
        original_dat["opensbliblock00"]["rho_B0"][halo_n:-halo_p, halo_n:-halo_p]
    )

    if args.halo_size == 0:
        rho_strided = numpy.array(
            strided_dat["opensbliblock00"]["rho_B0_strided"][:, :]
        )
    else:
        halo_n_new = int(args.halo_size)
        halo_p_new = int(args.halo_size)
        rho_strided = numpy.array(
            strided_dat["opensbliblock00"]["rho_B0_strided"][
                halo_n_new:-halo_p_new, halo_n_new:-halo_p_new
            ]
        )

    stride_i = int(args.stride_i)
    stride_j = int(args.stride_j)
    rho_original_strided = rho_original[::stride_i, ::stride_j]
    if rho_original_strided.shape != rho_strided.shape:
        print(
            "Shape mismatch original/strided = ",
            rho_original_strided.shape,
            rho_strided.shape,
        )
        return
    rho_diff = rho_original_strided - rho_strided

    fig, ax = pyplot.subplots(1, 3, figsize=(12, 5))
    ax[0].imshow(rho_original)
    ax[1].imshow(rho_strided)
    im = ax[2].imshow(rho_diff[:, :])
    fig.colorbar(im, ax=ax[2])
    ax[0].set_title("Original")
    ax[1].set_title("Strided")
    ax[2].set_title("Strided absolute difference")
    fig.tight_layout()
    pyplot.show()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("original_data")
    ap.add_argument("new_data")
    ap.add_argument("stride_i", default=1, type=int)
    ap.add_argument("stride_j", default=1, type=int)
    ap.add_argument("--halo_size", default=0, type=int)
    args = ap.parse_args()
    main(args)
