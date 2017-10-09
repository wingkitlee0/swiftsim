import numpy as np
import write_gadget as wg
import h5py as h5



class Particles(object):
    """
    Holder class for particle properties.
    """
    def __init__(self, meta):
        self.gravitymass = meta["gravitymass"]
        self.nparts = meta["nparts"]**2
        self.particlemass = meta["particlemass"]
        self.softening = meta["softening"]
        self.max = meta["max"]
        self.min = meta["min"]

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
        calculates keplerian orbit velocities for the, well, keplerian ring.
        """
        force_modifier = np.sqrt(self.gravitymass / (self.radii**2 + self.softening**2)**(3/2)) * self.radii 
        
        v_x = force_modifier * (-np.sin(self.theta))
        v_y = force_modifier * (np.cos(self.theta))
                    
        self.velocities = np.array([v_x, v_y, np.zeros_like(v_x)]).T

        return self.velocities


    def generate_ids(self):
        self.ids = np.arange(self.nparts)

        return self.ids


    def convert_polar_to_cartesian(self):
        x = self.radii * np.cos(self.theta)
        y = self.radii * np.sin(self.theta)

        self.positions = np.array([x, y, np.zeros_like(x)]).T

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

        try:
            self.smoothing = np.ma.array(self.smoothing, mask=mask).compressed()
            self.internalenergy = np.ma.array(self.internalenergy, mask=mask).compressed()
        except np.ma.core.MaskError:
            # QSP Fix has modified us
            self.smoothing = self.smoothing[:self.nparts]
            self.internalenergy = self.internalenergy[:self.nparts]

        try:
            self.masses = np.ma.array(self.masses, mask=mask).compressed()
        except np.ma.core.MaskError:
            # Not filled yet.
            pass

        self.radii = np.ma.array(self.radii, mask=mask).compressed()
        self.theta = np.ma.array(self.theta, mask=mask).compressed()

        self.nparts = len(self.radii)

        return


    def save_to_gadget(self, filename):
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

        ---------------------------------------------------------------------------
        """
        with h5.File(filename, "w") as handle:
            wg.write_header(
                handle,
                boxsize=8.,
                flag_entropy=0,
                np_total=np.array([self.nparts, 0, 0, 0, 0, 0]),
                np_total_hw=np.array([0, 0, 0, 0, 0, 0]),
                other={"MassTable" : np.array([self.particlemass, 0, 0, 0, 0, 0])}
            )

            wg.write_runtime_pars(
                handle,
                periodic_boundary=1,
            )

            wg.write_units(
                handle,
                current=1.,
                length=3.0856776e21,
                mass=1.9885e33,
                temperature=1.,
                time=3.085678e16,
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


sigma = np.vectorize(__sigma)


def generate_theta(r, theta_initial=0.):
    """
    Generate the theta associated with the particles based on their radii.

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


def gen_particles_grid(meta, range=(1, 7), centre_of_ring=(4, 4)):
    particles = Particles(meta)

    # Because we are using a uniform grid we actually use the same x and y
    # range for the initial particle setup.
    step = (range[1] - range[0])/meta["nparts"]
    
    x_values = np.arange(range[0], range[1], step)

    # These are 2d arrays which isn't actually that helpful.
    x, y = np.meshgrid(x_values, x_values)
    x = x.flatten()
    y = y.flatten()
    z = np.zeros_like(x)

    particles.positions = np.array([x, y, z]).T
    particles.radii = np.sqrt((x - centre_of_ring[0])**2 + (y - centre_of_ring[1])**2)
    particles.theta = np.arctan2(y - centre_of_ring[1], x - centre_of_ring[0])
    particles.exclude_particles((particles.softening, 100.))

    particles.densities = sigma(particles.radii)
    particles.calculate_velocities()
    particles.calculate_masses()

    particles.generate_ids()

    return particles


def gen_particles_spiral(meta, max_r=5., centre_of_ring=(4, 4)):
    particles = Particles(meta)

    # We follow the method of TODO
    m = (np.arange(particles.nparts) + 0.5)/particles.nparts
    r = max_r * m
    theta = generate_theta(r)

    particles.radii, particles.theta = QSP_fix(r, theta)
    particles.convert_polar_to_cartesian()
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
                    Hello
                    """
    )

    PARSER.add_argument(
        "-m",
        "--gravitymass",
        required=False,
        default=1.,
    )

    PARSER.add_argument(
        "-f",
        "--filename",
        required=False,
        default="initial_conditions.hdf5",
    )

    PARSER.add_argument(
        "-n",
        "--nparts",
        required=False,
        default=128,
    )

    PARSER.add_argument(
        "-p",
        "--particlemass",
        required=False,
        default=10.,
    )

    PARSER.add_argument(
        "-e",
        "--softening",
        required=False,
        default=0.05,
    )

    PARSER.add_argument(
        "-s",
        "--smoothing",
        required=False,
        default=0.89,
    )
    
    PARSER.add_argument(
        "-i",
        "--internalenergy",
        required=False,
        default=0.015
    )

    PARSER.add_argument(
        "-g",
        "--generationmethod",
        required=False,
        default="grid",
    )

    PARSER.add_argument(
        "-t",
        "--max",
        required=False,
        default=7.,
    )

    PARSER.add_argument(
        "-b",
        "--min",
        required=False,
        default=1.,
    )

    ### --- ### --- Argument Parsing --- ### --- ###

    ARGS = vars(PARSER.parse_args())

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
        "max": float(ARGS["max"]),
        "min": float(ARGS["min"]),
    }

    PARTICLES = gen_particles(META)

    PARTICLES.save_to_gadget(filename=ARGS["filename"])
