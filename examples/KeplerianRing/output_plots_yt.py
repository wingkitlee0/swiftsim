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
###############################################################################
"""

# Plotting script for the Keplerian Ring example.
# We use yt for the projection rendering of the ring,
# and then our own density as a function of radius calculation.

import matplotlib
matplotlib.use("Agg")

import yt
import h5py

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

from tqdm import tqdm
from make_movie import get_metadata


def get_axes_grid(figure):
    """
    Grab our axes grid.
    """
    gs = gridspec.GridSpec(2, 30)

    grid = []
    
    grid.append(figure.add_subplot(gs[0, 0:10]))
    grid.append(figure.add_subplot(gs[0, 10:20]))
    grid.append(figure.add_subplot(gs[0, 20:]))
    grid.append(figure.add_subplot(gs[1, 0:9]))
    grid.append(figure.add_subplot(gs[1, 11:20]))
    grid.append(figure.add_subplot(gs[1, 21:]))

    return grid


def get_yt_actual_data(plot):
    """
    Extracts the image data and colourmap from a yt plot.
    
    This is used to put on our own grid.
    """

    data = plot.plots["density"].image.get_array()
    cmap = plot.plots["density"].image.cmap

    return data, cmap


def chi_square(observed, expected):
    """
    The chi squared statistic.
    """

    return sum(((observed - expected)**2)/expected**2)


def load_data(filename):
    """
    Loads the data and extracts the relevant information for
    calculating the chi squared statistic and density(r) profiles.
    """

    with h5py.File(filename, "r") as file:
        boxsize = np.array(file["Header"].attrs["BoxSize"])

        # Check if z = 0 for all particles. If so we need to set the cetnre
        # to have z = 0 also.
        if np.sum(file["PartType0"]["Coordinates"][:, 2]) == 0:
            centre = [boxsize[0] / 2., boxsize[0] / 2., 0.]
        else:
            centre = boxsize / 2.

        radii = np.sqrt(np.sum(((file["PartType0"]["Coordinates"][...] - centre).T)**2, 0))
        masses = file["PartType0"]["Masses"][...]

    return radii, masses


def bin_density_r(radii, density, binrange, binnumber):
    """
    Bins the density as a funciton of radius.
    """

    bins = np.linspace(*binrange, binnumber)
    indicies = np.digitize(radii, bins)

    binned_masses = np.zeros(len(bins) - 1)
    
    for index, bin in enumerate(indicies):
        if bin >= len(bins) - 1:
            continue

        binned_masses[bin] += density[index]

    areas = [np.pi * (a**2 - b**2) for a, b in zip(bins[1:], bins[:-1])]
    binned_densities = binned_masses/areas
    
    return bins, binned_densities


def get_density_r(snapshot, filename="out", binrange=(0, 5), binnumber=50):
    """
    Gets the binned density as a function of radius.
    """
    snap = "{:04d}".format(snapshot)
    filename = f"{filename}_{snap}.hdf5"

    data = load_data(filename)

    return bin_density_r(*data, binrange, binnumber)


def get_derived_data(minsnap, maxsnap, filename="out"):
    """
    Gets the derived data from our snapshots, i.e. the
    density(r) profile and the chi squared (based on the
    difference between the minsnap and the current snapshot).
    """

    initial = get_density_r(minsnap, filename)
    densities = [get_density_r(snap)[1] for snap in tqdm(range(minsnap+1, maxsnap+1), desc="Densities")]
    densities = [initial[1]] + densities
    chisq = [chi_square(dens[5:20], initial[1][5:20]) for dens in tqdm(densities, desc="Chi Squared")]

    return initial[0], densities, chisq


def plot_chisq(ax, minsnap, maxsnap, chisq):
    """
    Plot the chisq(snapshot).
    """
    snapshots = np.arange(minsnap, maxsnap + 1)
    ax.plot(snapshots, np.array(chisq)/max(chisq))

    ax.set_xlabel("Snapshot Number")
    ax.set_ylabel("$\chi^2 / \chi^2_{{max}}$ = {:3.5f}".format(max(chisq)))

    return


def plot_density_r(ax, bins, densities, snaplist):
    """
    Make the density(r) plots.

    Densities is the _full_ list of density profiles, and
    snaplist is the ones that you wish to plot.
    """
    radii = [(x + y)/2 for x, y in zip(bins[1:], bins[:-1])]

    for snap in snaplist:
        index = snap - snaplist[0]
        ax.plot(radii, densities[index], label="Snapshot {:04d}".format(snap))

    ax.legend()
    ax.set_xlabel("Radius")
    ax.set_ylabel("Density")

    return

def plot_extra_info(ax, filename):
    """
    Plots all of the extra information on the final axis.

    Give it a filename of any of the snapshots.
    """

    metadata = get_metadata(filename)

    text = """
    Extra Information

    Box Size: {}

    Code
    
    Git Branch: {}
    Git Revision: {}
    Compiler: {} {}

    Hydro

    Scheme: {}
    Kernel $\eta$: {:4.4f}
    """.format(
        metadata['header']['BoxSize'][0]-1,
        metadata['code']['Git Branch'].decode("utf-8"),
        metadata['code']['Git Revision'].decode("utf-8"),
        metadata['code']['Compiler Name'].decode("utf-8"),
        metadata['code']['Compiler Version'].decode("utf-8"),
        metadata['hydro']['Scheme'].decode("utf-8"),
        metadata['hydro']['Kernel eta'][0],
    )
    
    ax.text(0.0, 0.1, text)
    ax.set_axis_off()

    return


def surface_density_plot(ax, snapnum, filename="out", density_limits=None):
    """
    Make the surface density plot (via yt).

    Also returns the max and minimum values for the density so these can
    be passed to the next call.
    """

    unit_base = {
        'length': (1.0, 'cm'),
        'velocity': (1.0, 'cm/s'),
        'mass': (1.0, 'g')
    }

    density_units = yt.units.gram / (yt.units.cm**3)
    
    snap = yt.load("{}_{:04d}.hdf5".format(filename, snapnum), unit_base=unit_base)

    max_density = snap.all_data()[("gas", "density")].max()
    min_density = snap.all_data()[("gas", "density")].min() + 1 * density_units
    
    new_density_limits = (min_density, max_density)

    if density_limits is None:
        density_limits = new_density_limits

    projection_plot = yt.ProjectionPlot(
        snap,
        "z",
        ("gas", "density"),
        width=6.
    )

    projection_plot.set_zlim("density", *density_limits)

    data = get_yt_actual_data(projection_plot)

    ax.imshow(data[0], cmap=data[1])
    ax.text(20, 50, "t = {:1.4f}".format(float(snap.current_time)), color='white')

    # We now want to remove all of the ticklabels.

    for axis in ['x', 'y']:
        ax.tick_params(
            axis=axis,          
            which='both',      
            bottom='off',      
            top='off',         
            left='off',
            right='off',
            labelleft='off',
            labelbottom='off'
        ) 

    return density_limits


if __name__ == "__main__":
    import sys

    filename = "out"
    snapshots = [int(x) for x in sys.argv[1:]]

    figure = plt.figure(figsize=(12, 10))
    axes = get_axes_grid(figure)
    figure.subplots_adjust(hspace=0, wspace=0)

    density_limits = None

    for snap, ax in zip(snapshots, axes[0:3]):
        density_limits = surface_density_plot(
            ax,
            snap,
            density_limits=density_limits
        )

    # Now we need to do the density(r) plot.

    # Derived data includes density profiles and chi squared
    derived_data = get_derived_data(snapshots[0], snapshots[2])

    plot_density_r(axes[3], derived_data[0], derived_data[1], snapshots)

    plot_chisq(axes[4], snapshots[0], snapshots[2], derived_data[2])

    plot_extra_info(axes[5], "out_0000.hdf5")

    plt.savefig("test.png")


