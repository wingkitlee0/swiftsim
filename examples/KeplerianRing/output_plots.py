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

from tqdm import tqdm


def load_data(filename, silent=True):
    if not silent:
        print(f"Loading data from {filename}")
    with h5.File(filename, "r") as file_handle:
        coords = file_handle['PartType0']['Coordinates'][...]
        time = float(file_handle['Header'].attrs['Time'])
        return coords, time


def rms(x):
    return np.sqrt(sum(x**2))


def get_metadata(filename):
    """ The metadata should be extracted from the first snapshot. """
    with h5.File(filename, "r") as file_handle:
        header = file_handle['Header']
        code = file_handle['Code']
        hydro = file_handle['HydroScheme']

        # we want to get the inner velocity of the ring.
        vel = file_handle['PartType0']['Velocities'][0]
        rad = file_handle['PartType0']['Coordinates'][0]
        period = 2 * np.pi * rms(rad-100) / rms(vel)

        return_values = {
            "header" : dict(header.attrs),
            "code" : dict(code.attrs),
            "period" : float(period),
            "hydro" : dict(hydro.attrs)
        }

    return return_values
        


def plot_single(number, scatter, text, metadata, ax, options=False):
    filename = "keplerian_ring_{:04d}.hdf5".format(number)
    coordinates, time = load_data(filename)

    text.set_text(
        "Time: {:1.2f} | Rotations {:1.2f}".format(
            time,
            time/metadata['period'],
        )
    )

    scatter.set_data(coordinates[:, 0], coordinates[:, 1])

    return scatter,


if __name__ == "__main__":
    import os

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    metadata = get_metadata("keplerian_ring_0000.hdf5")
    n_particle = metadata['header']['NumPart_Total'][0]

    # Initial plot setup

    scatter, = ax.plot([0]*n_particle, [0]*n_particle, ms=0.5, marker="o", linestyle="")
    ax.set_xlim(80, 120)
    ax.set_ylim(80, 120)

    time_text = ax.text(81, 81, "Time: {:1.2f} | Rotations {:1.2f}".format(
        0,
        0/metadata['period'],
    ))

    ax.text(81, 116, "Code: {} {} | {} {} \nHydro {}\n$\eta$={:1.4f}".format(
        metadata['code']['Git Branch'].decode("utf-8"),
        metadata['code']['Git Revision'].decode("utf-8"),
        metadata['code']['Compiler Name'].decode("utf-8"),
        metadata['code']['Compiler Version'].decode("utf-8"),
        metadata['hydro']['Scheme'].decode("utf-8"),
        metadata['hydro']['Kernel eta'][0],
    ))
    ax.set_title("Keplerian Ring Test")
    ax.set_xlabel("$x$ position")
    ax.set_ylabel("$y$ position")

    # Look for the number of files in the directory.
    i = 0
    while True:
        if os.path.isfile("keplerian_ring_{:04d}.hdf5".format(i)):
            i += 1
        else:
            break

        if i > 10000:
            break

    
    anim = anim.FuncAnimation(
        fig,
        plot_single,
        tqdm(np.arange(i)),
        fargs = [
            scatter,
            time_text,
            metadata,
            ax,
        ],
        interval=50,
        repeat_delay=3000,
        blit=True,
    )

    anim.save("keplerian_ring.mp4", dpi=int(640/8))
