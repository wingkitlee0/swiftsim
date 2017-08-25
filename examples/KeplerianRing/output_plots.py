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
        coords = file_handle['PartType0']['Coordinates'][...]
        time = float(file_handle['Header'].attrs['Time'])
        return coords, time

def rms(x):
    return sum(x**2)

def get_metadata(filename):
    """ The metadata should be extracted from the first snapshot. """
    with h5.File(filename, "r") as file_handle:
        header = file_handle['Header']
        code = file_handle['Code']
        hydro = file_handle['HydroScheme']

        # we want to get the inner velocity of the ring.
        vel = file_handle['PartType0']['Velocities'][0]
        rad = file_handle['PartType0']['Coordinates'][0]
        period = 2 * np.pi * rms(rad) / rms(vel)
        
        return_values = {
            "header" : dict(header.attrs),
            "code" : dict(code.attrs),
            "period" : float(period),
            "hydro" : dict(hydro.attrs)
        }

    return return_values
        


def plot_single(number, metadata, options=False):
    filename = "keplerian_ring_{:04d}.hdf5".format(number)
    coordinates, time = load_data(filename)

    plt.text(81, 81, "Time: {:1.2f} | Rotations {:1.2f}".format(
        time,
        time/metadata['period'],
    ))
    plt.text(81, 116, "Code: {} {} | {} {} \nHydro {}\n$\eta$={:1.4f}".format(
        metadata['code']['Git Branch'].decode("utf-8"),
        metadata['code']['Git Revision'].decode("utf-8"),
        metadata['code']['Compiler Name'].decode("utf-8"),
        metadata['code']['Compiler Version'].decode("utf-8"),
        metadata['hydro']['Scheme'].decode("utf-8"),
        metadata['hydro']['Kernel eta'][0],
    ))
    plt.title("Keplerian Ring Test")
    plt.xlabel("$x$ position")
    plt.ylabel("$y$ position")


    if options:
        return plt.scatter(coordinates[:, 0], coordinates[:, 1] ,**options)
    else:
        return plt.scatter(coordinates[:, 0], coordinates[:, 1])



if __name__ == "__main__":
    my_plots = []
    fig = plt.figure(figsize=(8, 8))
    metadata = get_metadata("keplerian_ring_0000.hdf5")

    i = 0
    while True:
        try:
            my_plots.append([plot_single(i, metadata, {"s" : 0.1, "c" : "b"})])
            plt.xlim(80, 120)
            plt.ylim(80, 120)
            i += 1
        except OSError:
            break

    anim = anim.ArtistAnimation(
        fig,
        my_plots,
        interval=50,
        repeat_delay=3000,
        blit=True,
    )

    anim.save("keplerian_ring.mp4", dpi=int(2048/8))
