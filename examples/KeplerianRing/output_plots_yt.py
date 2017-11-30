import yt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

def get_axes_grid(figure):
    """
    Grab our axes grid.
    """
    gs = gridspec.GridSpec(2, 3)

    grid = []
    
    grid.append(figure.add_subplot(gs[0, 0]))
    grid.append(figure.add_subplot(gs[0, 1]))
    grid.append(figure.add_subplot(gs[0, 2]))
    grid.append(figure.add_subplot(gs[1, 0]))
    grid.append(figure.add_subplot(gs[1, 1]))
    grid.append(figure.add_subplot(gs[1, 2]))

    return grid


def get_yt_actual_data(plot):
    data = plot.plots["density"].image.get_array()
    cmap = plot.plots["density"].image.cmap

    return data, cmap


def chi_square(observed, expected):
    return sum(((observed - expected)**2)/expected)

def get_density_r(snapshot, key="cell_mass"):
        snap = "{:04d}".format(snapshot)
        data[snap] = yt.load(f"{filename}_{snap}.hdf5", unit_base=unit_base, bounding_box=bbox)
        disks[snap] = data[snap].disk([4, 4, 0], [0, 0, 1], 3, 1) 
        prof = yt.Profile1D(disks[snap], "radius", 64, 0, 3, False, weight_field="cell_mass")
        prof.add_fields([key])

        density = prof[key][1:-1]

        return density


if __name__ == "__main__":
    import sys

    filename = "out"
    snapshots = ["{:04d}".format(int(x)) for x in sys.argv[1:]]

    figure = plt.figure(figsize=(20, 16))
    axes = get_axes_grid(figure)
    figure.subplots_adjust(hspace=0, wspace=0)

    unit_base = {
        'length': (1.0, 'cm'),
        'velocity': (1.0, 'cm/s'),
        'mass': (1.0, 'g')
    }

    bbox = [
        [-1, 9],
        [-1, 9],
        [-1, 9],
    ]

    density_units = yt.units.gram / (yt.units.cm**3)
    
    first_snap = yt.load(f"{filename}_{snapshots[0]}.hdf5", unit_base=unit_base)
    max_density = first_snap.all_data()[("gas", "density")].max()
    min_density = first_snap.all_data()[("gas", "density")].min() + 1 * density_units

    data = {}
    disks = {}

    for snap in snapshots:
        data[snap] = yt.load(f"{filename}_{snap}.hdf5", unit_base=unit_base, bounding_box=bbox)
        disks[snap] = data[snap].disk([4, 4, 0], [0, 0, 1], 3, 1) 

    for snap, ax in zip(snapshots, axes[0:3]):
        plot = yt.ProjectionPlot(
            data[snap],
            "z",
            ("gas", "density"),
            width=6.
        )

        plot.set_zlim("density", min_density, max_density)

        array_data, cmap = get_yt_actual_data(plot)
        ax.imshow(array_data, cmap=cmap, vmin=min_density, vmax=max_density)

    # Now we need to do the density(r) plot.

    density_r_data = {}
    r_data = {}
    key = "cell_mass"
    for snap in snapshots:
        plot = yt.ProfilePlot(disks[snap], "radius", key, "cell_mass", x_log=False, n_bins=64)
        x = plot.profiles[0].x
        density = plot.profiles[0][key]

        density_r_data[snap] = density
        r_data[snap] = x

    for snap in snapshots:
        axes[3].plot(r_data[snap], density_r_data[snap], label=f"Snapshot {snap}")

    # I guess we'll just have to hope they are all the same.
    # Which they are.
    axes[3].set_xlim(min(r_data[snapshots[0]]), max(r_data[snapshots[0]]))
    axes[3].set_xlabel(f"Distance from centre of ring (cm)")
    axes[3].set_ylabel(f"Cell mass (g)")
    axes[3].legend()

    zeroth_density = get_density_r(0)
    chisquare = []
    n_to_plot = 130
    for snapshot in range(120, n_to_plot):
        density = get_density_r(snapshot)

        chisquare.append(chi_square(density, zeroth_density))

    print(chisquare)
    snaps = range(120, n_to_plot)
    axes[4].plot(snaps, chisquare)
    axes[4].set_xlim(120, n_to_plot)
    axes[4].set_ylim(min(chisquare)*0.9, max(chisquare)*1.1)

    plt.savefig("test.png")
