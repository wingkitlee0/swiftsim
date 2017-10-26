import yt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

def get_axes_grid(figure):
    """
    Grab our axes grid.
    """
    grid = AxesGrid(figure, (0.075,0.075,0.85,0.85),
                nrows_ncols = (2, 3),
                axes_pad = 1.0,
                label_mode = "1",
                share_all = True,
                cbar_location="right",
                cbar_mode="single",
                cbar_size="3%",
                cbar_pad="0%")

    return grid


if __name__ == "__main__":
    import sys

    filename = "out"
    snapshots = ["{:04d}".format(int(x)) for x in sys.argv[1:]]

    figure = plt.figure(figsize=(20, 30))
    axes = get_axes_grid(figure)

    for snap, ax, cbax in zip(snapshots, axes[0:3], axes.cbar_axes[0:3]):
        data = yt.load(f"{filename}_{snap}.hdf5")

        plot = yt.ProjectionPlot(
            data,
            "z",
            ("gas", "density"),
            width=6.
        )

        p = plot.plots['density']
        p.figure = figure
        p.axes = ax
        p.cax = cbax

        plot._setup_plots()

    plt.savefig("test.png")
