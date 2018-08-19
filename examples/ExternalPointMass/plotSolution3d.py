
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import h5py
import argparse

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('infile',type=str,help="hdf5 file to plot")
    args = parser.parse_args()


    sim = h5py.File(args.infile, "r")
    x = sim["PartType1/Coordinates"][:,0]
    y = sim["PartType1/Coordinates"][:,1]
    z = sim["PartType1/Coordinates"][:,2]
    vx = sim["PartType1/Velocities"][:,0]
    vy = sim["PartType1/Velocities"][:,1]
    vz = sim["PartType1/Velocities"][:,2]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x,y,z)
    plt.show()




