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


import matplotlib.pyplot as plt
from generate_ics import generate_particles


if __name__ == "__main__":
    # Check that keplerian velocities look correct.
    x, y, vx, vy = generate_particles(100, 10, 2.5, 1000)

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)

    ax.quiver(x, y, vx, vy)
    ax.set_xlim(-20, 20)
    ax.set_ylim(-20, 20)

    fig.show()
    print("Press enter to quit.")
    input()  # keep the figure alive
