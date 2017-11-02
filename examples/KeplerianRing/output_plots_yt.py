import yt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import AxesGrid

def get_axes_grid(figure):
    """
    Grab our axes grid.
    """
    grid = AxesGrid(figure, (0.075,0.075,0.85,0.85),
                nrows_ncols = (2, 3),
                axes_pad = 0, 
                label_mode = "1",
                share_all = True,
                cbar_location="right",
                cbar_mode="edge",
                cbar_size="3%",
                cbar_pad="0%",
                aspect=False,)

    return grid


if __name__ == "__main__":
    import sys

    filename = "out"
    snapshots = ["{:04d}".format(int(x)) for x in sys.argv[1:]]

    figure = plt.figure(figsize=(4, 4))
    axes = get_axes_grid(figure)
    figure.subplots_adjust(hspace=0, wspace=0)

    unit_base = {
        'length': (1.0, 'cm'),
        'velocity': (1.0, 'cm/s'),
        'mass': (1.0, 'g')
    }

    density_units = yt.units.gram / (yt.units.cm**3)
    
    first_snap = yt.load(f"{filename}_{snapshots[0]}.hdf5", unit_base=unit_base)
    max_density = first_snap.all_data()[("gas", "density")].max()
    min_density = first_snap.all_data()[("gas", "density")].min() + 1 * density_units

    density_cbar = axes.cbar_axes[0]

    data = {}
    disks = {}

    for snap in snapshots:
        data[snap] = yt.load(f"{filename}_{snap}.hdf5", unit_base=unit_base)
        disks[snap] = data[snap].disk([8, 8, 0], [0, 0, 1], 4, 1) 

    for snap, ax in zip(snapshots, axes[0:3]):
        plot = yt.ProjectionPlot(
            data[snap],
            "z",
            ("gas", "density"),
            width=6.
        )

        plot.set_zlim("density", min_density, max_density)

        p = plot.plots['density']
        p.figure = figure
        p.axes = ax
        p.cax = density_cbar
        ax.set_xlim(-3, 3)
        ax.set_ylim(-3, 3)

        plot._setup_plots()

    # Now we need to do the density(r) plot.

    density_r_data = {}
    r_data = {}
    for snap in snapshots:
        plot = yt.ProfilePlot(disks[snap], "particle_radius", "density", "particle_mass")
        x = plot.profiles[0].x
        density = plot.profiles[0]["density"]

        density_r_data[snap] = density
        r_data[snap] = x

        x_units = plot.profiles[0]

    for snap in snapshots:
        axes[3].plot(r_data[snap], density_r_data[snap], label=snap)

    # I guess we'll just have to hope they are all the same.
    # Which they are.
    axes[3].set_xlim(max(density_r_data[snapshots[0]]), min(density_r_data[snapshots[0]]))
    axes[3].set_xlabel(f"Density ${r_data[snapshots[0]].units}$")
    axes[3].set_ylabel(f"Displacement from ring centre ${density_r_data[snapshots[0]].units}$")
    axes[3].legend()

    plt.savefig("test.png")
