
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import h5py
import argparse

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('infile',type=str,help="hdf5 file to plot")
    parser.add_argument('-pt',type=int,default=0,help="particle type (int)")
    parser.add_argument('-twod',action='store_true',help="2d plot of particles")
    args = parser.parse_args()

    coordname = "PartType%i/Coordinates" % args.pt
    velname = "PartType%i/Velocities" % args.pt

    sim = h5py.File(args.infile, "r")


    x = sim[coordname][:,0]
    y = sim[coordname][:,1]
    z = sim[coordname][:,2]
    vx = sim[velname][:,0]
    vy = sim[velname][:,1]
    vz = sim[velname][:,2]

    if args.twod:
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect="equal")
        ax.scatter(x,y,alpha=0.1,s=10)
    else:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x,y,z)
    plt.show()




