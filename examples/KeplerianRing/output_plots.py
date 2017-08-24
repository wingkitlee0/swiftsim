"""
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2017
#
# Josh Borrow (joshua.borrow@durham.ac.uk)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
# -----------------------------------------------------------------------------
#
# This program creates test plots for the initial condition generator provided
# for the Keplerian Ring example.
#
###############################################################################
"""


import matplotlib.pyplot as plt
import matplotlib.animation as anim
import numpy as np
import h5py as h5


def load_data(filename):
    with h5.File(filename, "r") as file_handle:
        return file_handle['PartType0']['Coordinates']


def plot_single(number, ax, options=False):
    filename = "keplerian_ring_{:04d}".format(number)
    coordinates = load_data(filename)

    if options:
        ax.scatter(coordinates[:, 0], coordinates[:, 1] ,**options)
    else:
        ax.scatter(coordinates[:, 0], coordinates[:, 1])

    return ax


if __name__ == "__main__":
    my_plots = []
    fig, ax = plt.subplots()

    for i in range(12):
        my_plots.append(plot_single(i, ax))

    ax.set_xlim(-100, 100)
    ax.set_ylim(-100, 100)

    anim = anim.ArtistAnimation(
        fig,
        my_plots,
        interval=50,
        repeat_delay=3000,
        blit=True,
    )

    anim.save("keplerian_ring.mp4")
