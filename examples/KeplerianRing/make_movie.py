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

import matplotlib
matplotlib.use("Agg")


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


def rotation_velocity_at_r(r, params):
    """
    Gets the rotation velocity at a given radius r by assuming it is keplerian.

    Assumes we are in cgs units, which may one day be our downfall.
    """

    unit_length = float(params[r"InternalUnitSystem:UnitLength_in_cgs"])

    if unit_length != 1.:
        print(f"Your unit length: {unit_length}")
        raise InternalUnitSystemError(
            "This function is only created to handle CGS units."
        )

    central_mass = float(params["PointMassPotential:mass"])
    G = 6.674e-8

    v = np.sqrt( G * central_mass / r)

    return v


def get_rotation_period_at_r(r, params):
    """
    Gets the rotation period at a given radius r, assuming a keplerian
    orbit.
    """
    v = rotation_velocity_at_r(r, params)

    return 2*np.pi / v


def get_metadata(filename, r=1):
    """ The metadata should be extracted from the first snapshot. """
    with h5.File(filename, "r") as file_handle:
        header = file_handle['Header'].attrs
        code = file_handle['Code'].attrs
        hydro = file_handle['HydroScheme'].attrs
        params = file_handle['Parameters'].attrs

        period = get_rotation_period_at_r(r, params)

        return_values = {
            "header" : dict(header),
            "code" : dict(code),
            "period" : float(period),
            "hydro" : dict(hydro),
            "params" : dict(params)
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
    ax.set_xlim(0, metadata['header']['BoxSize'][0])
    ax.set_ylim(0, metadata['header']['BoxSize'][1])

    offset = 0.25
    time_text = ax.text(offset, offset, "Time: {:1.2f} | Rotations {:1.2f}".format(
        0,
        0/metadata['period'],
    ))

    ax.text(offset, metadata['header']['BoxSize'][0]-offset-0.5, "Code: {} {} | {} {} \nHydro {}\n$\eta$={:1.4f}".format(
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
