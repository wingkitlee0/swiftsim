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
# This program is a set of helper functions that write particle data to a
# GADGET .hdf5 file *handle*.
#
###############################################################################
"""

def write_block(f, part_type, pos, vel, id, m=False, U=False, rho=False, hsml=False, pot=False, acc=False, dAdt=False, dt=False):
    """
    Write a block of data (i.e. a set of particle types).

    @param: f | file handle
        - the hdf5 file handle to write to

    @param: part_type | int
        - the GADGET particle type number. See the GADGET documentation.

    @param: pos | array-like
        - the array of particle positions with shape (n_particles, 3)

    @param: vel | array-like
        - the array of particle velocities with shape (n_particles, 3)
    
    @param: id | array-like
        - array of particle ids with shape (n_particles, 3)

    @param: m | array-like | optional
        - array of particle masses (only useful if particle masses change with
          time, otherwise just set the mass in the header).

    @param: U | array-like | optional
        - array of the internal energies associated with the particles.

    @param: rho | array-like | optional
        - array of the densities of the particles.

    @param: hsml | array-like | optional
        - array of the smoothing lengths of the particles (again, this can be
          set for all particles elsewhere).
    
    ...
    """
    part_type = f.create_group("PartType{}".format(part_type))
    n_particles = len(id)

    # Create datasets
    positions = type.create_dataset("Coordinates", (n_particles, 3), data=pos)
    velocities = type.create_dataset("Velocities", (n_particles, 3))
    ids = type.create_dataset("ParticleIDs", (n_particles,))

    # Oh gosh this is very wasteful code but it *works*.
    if m: #.any()
        masses = type.create_dataset("Masses", (n_particles, ))
    if U: #.any()
        internal_energies = type.create_dataset("InternalEnergy", (n_particles, ))
    if rho: #.any()
        density = type.create_dataset("Density", (n_particles, ))
    if hsml: #.any()
        smoothing_length = type.create_dataset("SmoothingLength", (n_particles, ))
    if pot: #.any()
        potential = type.create_dataset("Potential", (n_particles, ))
    if acc: #.any()
        accelerations = type.create_dataset("Acceleration", (n_particles, 3))
    if dAdt: #.any()
        rate_of_change_of_entropy = type.create_dataset("RateOfChangeOfEntropy", (n_particles, ))
    if dt: #.any()
        timestep = type.create_dataset("TimeStep", (n_particles, ))

    # Now fill

    velocities[...] = vel
    ids[...] = id

    if m: #.any()
        masses[...] = m
    if U: #.any()
        internal_energies[...] = U
    if rho: #.any()
        density[...] = rho
    if hsml: #.any()
        smoothing_length[...] = hsml
    if pot: #.any()
        potential[...] = pot
    if acc: #.any()
        accelerations[...] = acc
    if dAdt: #.any()
        rate_of_change_of_entropy[...] = dAdt
    if dt: #.any()
        timestep[...] = dt

    return f


def write_head(f, npart, massarr, time, z=False, flagsfr=False, flagfeedback=False, nall=False, flagcooling=False, numfiles=1, omega0=0.3, omegalambda=0.7, hubbleparam=0.7, flagage=False, flagmetals=False, nallhw=0, flagentrics=False):
    header = f.create_group('Header')

    header.attrs['NumPart_ThisFile'] = npart
    header.attrs['MassTable'] = massarr
    header.attrs['Time'] = time
    header.attrs['NumFilesPerSnapshot'] = numfiles
    header.attrs['Omega0'] = omega0
    header.attrs['OmegaLambda'] = omegalambda
    header.attrs['HubbleParam'] = hubbleparam
    header.attrs['NumPart_Total_HW'] = nallhw

    if z: #.any()
        header.attrs['Redshift'] = z
    if flagsfr: #.any()
        header.attrs['Flag_Sfr'] = flagsfr
    if flagfeedback: #.any()
        header.attrs['Flag_Feedback'] = flagfeedback
    if nall: #.any()
        header.attrs['NumPart_Total'] = nall
    else:
        header.attrs['NumPart_Total'] = npart
    if flagcooling: #.any()
        header.attrs['Flag_Cooling'] = flagcooling
    if flagage: #.any()
        header.attrs['Flag_StellarAge'] = flagage
    if flagmetals: #.any()
        header.attrs['Flag_Metals'] = flagmetals
    if flagentrics: #.any()
        header.attrs['Flag_Entropy_ICs'] = flagentrics

    return f