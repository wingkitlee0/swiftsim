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


import numpy as np
import write_gadget as wg
import h5py as h5


class Particles(object):
    """
    Holder class for particle properties. Also contains some methods to e.g.
    set their keplerian velocities based on their positions. These properties
    are set using the 'generationmethod' functions below.
    """
    def __init__(self, meta):
        self.gravitymass = meta["gravitymass"]
        self.nparts = meta["nparts"]**2
        self.particlemass = meta["particlemass"]
        self.softening = meta["softening"]
        self.boxsize = meta["boxsize"]

        self.smoothing = np.zeros(self.nparts) + meta["smoothing"]
        self.internalenergy = np.zeros(self.nparts) + meta["internalenergy"]

        self.positions = np.array([])
        self.radii = np.array([])
        self.theta = np.array([])
        self.velocities = np.array([])
        self.ids = np.array([])
        self.densities = np.array([])
        self.masses = np.array([])

        return


    def calculate_masses(self):
        """
        Calculate the individual masses for the particles such that we have
        a uniform thickness disk, given their densities.
        """
        mass_factor = self.particlemass / self.densities.max()
        self.masses = self.densities * mass_factor

        return self.masses
   

    def calculate_velocities(self):
        """
        Calculates keplerian orbit velocities for the, well, keplerian ring.
        Requires that radius and theta are set already.
        """
        force_modifier = np.sqrt(self.gravitymass / (self.radii**2 + self.softening**2)**(3/2)) * self.radii 
        
        v_x = force_modifier * (-np.sin(self.theta))
        v_y = force_modifier * (np.cos(self.theta))
                    
        self.velocities = np.array([v_x, v_y, np.zeros_like(v_x)]).T

        return self.velocities


    def generate_ids(self):
        """
        Generate consecutive IDs from 0 based on the number of particles
        currently in the object.
        """
        self.ids = np.arange(self.nparts)

        return self.ids


    def convert_polar_to_cartesian(self, centre_of_ring=(5, 5), boxsize):
        """
        Converts self.radii, self.theta to self.positions.
        """
        x = self.radii * np.cos(self.theta) + centre_of_ring[0]
        y = self.radii * np.sin(self.theta) + centre_of_ring[1]
        z = np.zeros_like(x) + boxsize/2

        self.positions = np.array([x, y, z]).T

        return self.positions


    def exclude_particles(self, range):
        """
        Exclude all particles that are *not* within range (of radii).
        """
        mask = np.logical_or(
            self.radii < range[0],
            self.radii > range[1],
        )

        x  = np.ma.array(self.positions[:, 0], mask=mask).compressed()
        y  = np.ma.array(self.positions[:, 1], mask=mask).compressed()
        z  = np.ma.array(self.positions[:, 2], mask=mask).compressed()

        self.positions = np.array([x, y, z]).T
        
        try:
            v_x  = np.ma.array(self.velocities[:, 0], mask=mask).compressed()
            v_y  = np.ma.array(self.velocities[:, 1], mask=mask).compressed()
            v_z  = np.ma.array(self.velocities[:, 2], mask=mask).compressed()

            self.velocities = np.array([v_x, v_y, v_z])
        except IndexError:
            # We haven't filled them yet anyway...
            pass
    
        try:
            self.ids = np.ma.array(self.ids, mask=mask).compressed()
        except np.ma.core.MaskError:
            # Not filled yet.
            pass

        try:
            self.densities = np.ma.array(self.densities, mask=mask).compressed()
        except np.ma.core.MaskError:
            # Not filled yet.
            pass

        # QSP Fix has modified us, so first we need to chop off extras.
        # Then, as they are all the same, we don't need to remove 'specific'
        # values, and so we can just chop off the ends.
        self.smoothing = self.smoothing[:len(x)]
        self.internalenergy = self.internalenergy[:len(x)]

        try:
            self.masses = np.ma.array(self.masses, mask=mask).compressed()
        except np.ma.core.MaskError:
            # Not filled yet.
            pass

        self.radii = np.ma.array(self.radii, mask=mask).compressed()
        self.theta = np.ma.array(self.theta, mask=mask).compressed()

        self.nparts = len(self.radii)

        return


    def save_to_gadget(self, filename, boxsize=10.):
        """ Save the particle data to a GADGET .hdf5 file.

        @param: filename | string
        - filename of the hdf5 file to save.

        @param: x_i | array-like
        - x positions of the particles

        @param: y_i | array-like
        - y positions of the particles

        @param: v_x_i | array-like
        - x velocities of the particles

        @param: v_y_i | array-like
        - y velocities of the particles

        @param: hsml | float
        - smoothing length of the particles.

        @param: pm | float
        - mass of the particles.

        -----------------------------------------------------------------------
        """
        with h5.File(filename, "w") as handle:
            wg.write_header(
                handle,
                boxsize=boxsize,
                flag_entropy=0,
                np_total=np.array([self.nparts, 0, 0, 0, 0, 0]),
                np_total_hw=np.array([0, 0, 0, 0, 0, 0]),
                other={
                    "MassTable" : np.array([self.particlemass, 0, 0, 0, 0, 0]),
                    "Time" : 0,
                }
            )

            wg.write_runtime_pars(
                handle,
                periodic_boundary=1,
            )

            wg.write_units(
                handle,
                current=1.,
                length=1.,
                mass=1,
                temperature=1.,
                time=1.,
            )

            wg.write_block(
                handle,
                0,  # gas
                self.positions,
                self.velocities,
                self.ids,
                mass=self.masses,
                int_energy=self.internalenergy,
                smoothing=self.smoothing,
                other={"Density": self.densities},
            )

        return


def __sigma(r):
    """
    Density distribution of the ring, this comes directly from Hopkins 2015.
    """
    if r < 0.5:
        return 0.01 + (r/0.5)**3
    elif r <= 2:
        return 1.01
    else:
        return 0.01 + (1 + (r-2)/0.1)**(-3)


# This is required because of the if, else statement above that does not
# play nicely with our numpy arrays.
sigma = np.vectorize(__sigma)


def generate_theta(r, theta_initial=0.):
    """
    Generate the theta associated with the particles based on their radii.
    This uses the method from The effect of Poisson noise on SPH calculations,
    Annabel Cartwright, Dimitris Stamatellos and Anthony P. Whitworth 2009.

    @param: r | float / array-like
        - the radii of the particles.

    @param: theta_initial | float | optional
        - initial radius to start the generation at.

    ---------------------------------------------------------------------------

    @return: theta_i | numpy.array
        - angles associated with the particles in the plane.
    """
    radii_fraction = r[:-1] / r[1:]
    d_theta_i = np.sqrt(2 * np.pi * (1 - radii_fraction))

    theta_i = np.empty_like(r)
    theta_i[0] = theta_initial

    # Need to do this awful loop because each entry relies on the one before it.
    # Unless someone knows how to do this in numpy.
    for i in range(len(d_theta_i)):  # first is theta_initial
        theta_i[i+1] = theta_i[i] + d_theta_i[i]

    return theta_i


def QSP_fix(r_i, theta_i):
    """
    The start and end of the disk will have the end of the spiral there. That
    is no good and introduces some shear forces, so we need to move them to
    concentric circles. We'll also add some extra particles on the final
    'layer' of this giant onion to ensure that all of the circles are complete.

    @param: r_i | numpy.array
        - the original r_i generated from inverse_gaussian (and perhaps
          afterwards masked).

    @param: theta_i | numpy.array
        - the original theta_i generated from generate_theta_i.

    ---------------------------------------------------------------------------

    @return r_i_fixed | numpy.array
        - the fixed, concentric circle-like r_i. Note that these arrays do not
          necessarily have the same length as the r_i, theta_i that are
          passed in and you will have to re-calculate the number of particles
          in the system.

    @return theta_i_fixed | numpy.array
        - the fixed, concentric circle-like theta_i
    """

    # Thankfully, the delta_thetas are not wrapped (i.e. they keep on going
    # from 2 pi -- we need to split the arrays into 'circles' by splitting
    # the theta_i array every 2 pi.

    rotations = 1
    circles = []

    these_r_i = []
    these_theta_i = []

    for radius, theta in zip(r_i, theta_i):
        if theta > rotations * 2 * np.pi:
            circles.append([these_r_i, these_theta_i])
            these_r_i = []
            these_theta_i = []
            rotations += 1

        these_r_i.append(radius)
        these_theta_i.append(theta)

    # Now we need to do the averaging.
    # We want to have all particles in each circle a fixed radius away from the
    # centre, as well as having even spacing between each particle. The final
    # ring may be a bit dodgy still, but we will see.
    
    r_i_fixed = []
    theta_i_fixed = []

    for circle in circles:
        n_particles = len(circle[0])
        radius = sum(circle[0]) / n_particles
        theta_initial = circle[1][0]

        theta_sep = 2 * np.pi / n_particles
        
        theta = [t * theta_sep for t in range(n_particles)]
        radii = [radius] * n_particles

        r_i_fixed += radii
        theta_i_fixed += theta

    return np.array(r_i_fixed), np.array(theta_i_fixed)


def gen_particles_grid(meta):
    """
    Generates particles on a grid and returns a filled Particles object.
    """
    particles = Particles(meta)
    range = (0, meta["boxsize"])
    centre_of_ring = (meta["boxsize"]/2., meta["boxsize"]/2.)

    # Because we are using a uniform grid we actually use the same x and y
    # range for the initial particle setup.
    step = (range[1] - range[0])/meta["nparts"]
    
    x_values = np.arange(0, range[1] - range[0], step)

    # These are 2d arrays which isn't actually that helpful.
    x, y = np.meshgrid(x_values, x_values)
    x = x.flatten() + centre_of_ring[0] - (range[1] - range[0])/2
    y = y.flatten() + centre_of_ring[1] - (range[1] - range[0])/2
    z = np.zeros_like(x) + meta["boxsize"]/2

    particles.positions = np.array([x, y, z]).T
    particles.radii = np.sqrt((x - centre_of_ring[0])**2 + (y - centre_of_ring[1])**2)
    particles.theta = np.arctan2(y - centre_of_ring[1], x - centre_of_ring[0])
    particles.exclude_particles((particles.softening, 100.))

    particles.densities = sigma(particles.radii)
    particles.calculate_velocities()
    particles.calculate_masses()

    particles.generate_ids()

    return particles


def gen_particles_spiral(meta):
    """
    Generates particles on concentric circles and returns a filled Particles
    object. Based on Cartwright, Stamatellos & Whitworth (2009).
    """
    particles = Particles(meta)
    centre_of_ring = (meta["boxsize"]/2., meta["boxsize"]/2.)
    max_r = meta["boxsize"]/2.


    m = (np.arange(particles.nparts) + 0.5)/particles.nparts
    r = max_r * m
    theta = generate_theta(r)

    particles.radii, particles.theta = QSP_fix(r, theta)
    particles.convert_polar_to_cartesian(centre_of_ring, meta["boxsize"])
    particles.nparts = len(particles.radii)
    
    particles.exclude_particles((particles.softening, 100.))
    
    particles.densities = sigma(particles.radii)
    particles.calculate_velocities()
    particles.calculate_masses()

    particles.generate_ids()

    return particles


if __name__ == "__main__":
    import argparse as ap

    # TODO: Add these descriptions
    PARSER = ap.ArgumentParser(
        description="""
                    Initial conditions generator for the Keplerian Ring
                    example. It has sensible defaults for GIZMO, but if you
                    wish to run the example with SPH you sould use
                    --generationmethod spiral.
                    """
    )

    PARSER.add_argument(
        "-m",
        "--gravitymass",
        help="""
             GM for the central point mass. Default: 1.
             """,
        required=False,
        default=1.,
    )

    PARSER.add_argument(
        "-f",
        "--filename",
        help="""
             Filename for your initial conditions.
             Default: initial_conditions.hdf5.
             """,
        required=False,
        default="initial_conditions.hdf5",
    )

    PARSER.add_argument(
        "-n",
        "--nparts",
        help="""
             Square-root of the number of particles, i.e. the default
             nparts=128 leads to a ring with 128^2 particles in it.
             Note that for the spiral generation method this number is 
             somewhat approximate.
             """,
        required=False,
        default=128,
    )

    PARSER.add_argument(
        "-p",
        "--particlemass",
        help="""
             Maximum mass of the gas particles. Default: 10.
             """,
        required=False,
        default=10.,
    )

    PARSER.add_argument(
        "-e",
        "--softening",
        help="""
             Gravitational softening for the centre of the ring. We also
             exclude all particles within this radius from the centre of the
             ring. Default: 0.05.
             """,
        required=False,
        default=0.05,
    )

    PARSER.add_argument(
        "-s",
        "--smoothing",
        help="""
             Initial smoothing length for all of the particles.
             Default: 0.89.
             """,
        required=False,
        default=0.89,
    )
    
    PARSER.add_argument(
        "-i",
        "--internalenergy",
        help="""
             Initial internal energy for all of the particles. Note that this
             should be low to ensure that the ring is very cold (see Inviscid
             SPH for details). Default: 0.015.
             """,
        required=False,
        default=0.015
    )

    PARSER.add_argument(
        "-g",
        "--generationmethod",
        help="""
             Generation method for the particles. Choose between grid and
             spiral, where spiral ensures that the particles are generated
             in a way that minimises the energy in SPH. For more details on
             this method see Cartwright, Stamatellos & Whitworth (2009).
             Default: grid.
             """,
        required=False,
        default="grid",
    )

    PARSER.add_argument(
        "-b",
        "--boxsize",
        help="""
             The box size. Default: 10.
             """,
        required=False,
        default=10.
    )

    PARSER.add_argument(
        "-a",
        "--angle",
        help="""
             Angle of incline on the disk (degrees). Default: 0.
             """,
        required=False,
        default=0.
    )


    ### --- ### --- Argument Parsing --- ### --- ###

    ARGS = vars(PARSER.parse_args())

    if ARGS["angle"]:
        print("Varing angles are not yet implemented.")
        raise NotImplementedError

    if ARGS["generationmethod"] == "grid":
        gen_particles = gen_particles_grid
    elif ARGS["generationmethod"] == "spiral":
        gen_particles = gen_particles_spiral
    else:
        print(
            "ERROR: {} is an invalid generation method. Exiting.".format(
                ARGS["generationmethod"]
            )
        )
        exit(1)

    META = {
        "gravitymass": float(ARGS["gravitymass"]),
        "nparts": int(ARGS["nparts"]),
        "particlemass": float(ARGS["particlemass"]),
        "smoothing": float(ARGS["smoothing"]),
        "softening": float(ARGS["softening"]),
        "internalenergy": float(ARGS["internalenergy"]),
        "boxsize": float(ARGS["boxsize"]),
        "angle" : float(ARGS["angle"])
    }

    PARTICLES = gen_particles(META)

    PARTICLES.save_to_gadget(
        filename=ARGS["filename"],
        boxsize=ARGS["boxsize"],
    )

