"""Compare reference data against output data for a wave equation OpenSBLI
simulation.
"""

import argparse
import sys

import h5py
import numpy
from matplotlib import pyplot
import plotly.graph_objects as go
from IPython import embed


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

    rho_old = numpy.array(original["opensbliblock00"]["rho_B0"][halo_n:-halo_p, halo_n:-halo_p, halo_n:-halo_p])
    rho_new = numpy.array(new["opensbliblock00"]["rho_B0"][halo_n:-halo_p, halo_n:-halo_p, halo_n:-halo_p])

    print(rho_old.shape, rho_new.shape)

    rho_old=rho_old[:,:,0]
    rho_new=rho_new[:,:,0]
    diff=rho_new-rho_old
    print(rho_old)
    print(rho_new)
    print(diff)

    ax.imshow(diff)

    # X, Y, Z = numpy.mgrid[0:rho.shape[0], 0:rho.shape[1], 0:rho.shape[2]]

    # fig = go.Figure(data=go.Volume(
    #     x=X.flatten(),
    #     y=Y.flatten(),
    #     z=Z.flatten(),
    #     value=rho.flatten(),
    #     opacity=0.8,
    #     surface_count=15,
    #     colorscale='viridis',
    #     caps=dict(x_show=False, y_show=False, z_show=False)  # Hide end-caps
    # ))

    pyplot.show()

    # ax.set_xlabel("$x$")
    # ax.set_ylabel(r"$\phi(x)$")
    # ax.legend(loc="upper right")
    # pyplot.show()

    return fig, ax

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("original_data")
    ap.add_argument("new_data")
    args = ap.parse_args()
    fig, ax = main(args)
