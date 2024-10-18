"""Compare reference data against output data for a wave equation OpenSBLI
simulation.
"""

import argparse
import sys

import h5py
from matplotlib import pyplot


def main(args):
    try:
        original = h5py.File(args.original_data, "r")
    except IOError:
        print(f"Unable to open: {args.original_data}")
        sys.exit(1)
    try:
        new = h5py.File(args.new_data, "r")
    except IOError:
        print(f"Unable to open: {args.new_data}")
        sys.exit(1)

    halo_n = 5
    halo_p = 5
    fig, ax = pyplot.subplots(1, 1, figsize=(8, 5))

    x = original["opensbliblock00"]["x0_B0"] [halo_n:-halo_p]
    y = original["opensbliblock00"]["phi_B0"][halo_n:-halo_p]
    ax.plot(x, y, label="Original")
    print(x.size, x)
    print(y.size, y)

    x = new["opensbliblock00"]["x0_B0"][:]
    y = new["opensbliblock00"]["phi_B0"][:]
    ax.plot(x, y, label="New (strided)")
    print(x.size, x)
    print(y.size, y)

    ax.set_xlabel("$x$")
    ax.set_ylabel(r"$\phi(x)$")
    ax.legend(loc="upper right")

    pyplot.show()

    return fig, ax

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("original_data")
    ap.add_argument("new_data")
    args = ap.parse_args()
    fig, ax = main(args)
