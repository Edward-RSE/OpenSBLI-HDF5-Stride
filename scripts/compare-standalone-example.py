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
    strided_name = Path(args.new_data).stem

    halo_n = 5
    halo_p = 5
    rho_original = numpy.array(
        original_dat["opensbliblock00"]["rho_B0"][
            halo_n:-halo_p, halo_n:-halo_p, halo_n:-halo_p
        ]
    )

    if args.halo_size == 0:
        rho_strided = numpy.array(strided_dat["opensbliblock00"]["rho_B0"][:, :, :])
    else:
        halo_n_new = int(args.halo_size)
        halo_p_new = int(args.halo_size)
        rho_strided = numpy.array(
            strided_dat["opensbliblock00"]["rho_B0"][
                halo_n_new:-halo_p_new, halo_n_new:-halo_p_new, halo_n_new:-halo_p_new
            ]
        )

    stride_i = int(args.stride_i)
    stride_j = int(args.stride_j)
    stride_k = int(args.stride_k)
    rho_original = rho_original[::stride_i, ::stride_j, ::stride_k]
    if rho_original.shape != rho_strided.shape:
        print(
            "Shape mismatch original/strided = ", rho_original.shape, rho_strided.shape
        )
        return
    rho_diff = rho_original - rho_strided

    for i in range(rho_original.shape[2]):
        fig, ax = pyplot.subplots(1, 1, figsize=(8, 5))
        diff = rho_diff[:, :, i]
        im = ax.imshow(diff)
        _cbar = fig.colorbar(im, ax=ax)
        ax.set_title("Absolute difference")
        path = Path(__file__).parent
        pyplot.savefig(f"{path}/fig/{strided_name}-taylor-green-" + str(i) + ".png")
        pyplot.close()
        if numpy.sum(diff) != 0:
            print(i, diff)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("original_data")
    ap.add_argument("new_data")
    ap.add_argument("stride_i", default=1, type=int)
    ap.add_argument("stride_j", default=1, type=int)
    ap.add_argument("stride_k", default=1, type=int)
    ap.add_argument("--halo_size", default=0, type=int)
    args = ap.parse_args()
    main(args)
