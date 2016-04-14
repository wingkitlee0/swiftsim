/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Tom Theuns (tom.theuns@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

/* Config parameters. */
#include "../config.h"

/* This object's header. */
#include "potentials.h"

/**
 * @brief Initialises the external potential properties in the internal system
 * of units.
 *
 * @param parameter_file The parsed parameter file
 * @param us The current internal system of units
 * @param potential The external potential properties to initialize
 */
void potential_init(const struct swift_params* parameter_file,
                    struct UnitSystem* us,
                    struct external_potential* potential) {

#ifdef EXTERNAL_POTENTIAL_POINTMASS

  potential->point_mass.x =
      parser_get_param_double(parameter_file, "PointMass:position_x");
  potential->point_mass.y =
      parser_get_param_double(parameter_file, "PointMass:position_y");
  potential->point_mass.z =
      parser_get_param_double(parameter_file, "PointMass:position_z");
  potential->point_mass.mass =
      parser_get_param_double(parameter_file, "PointMass:mass");

#endif /* EXTERNAL_POTENTIAL_POINTMASS */
}

/**
 * @brief Prints the properties of the external potential to stdout.
 *
 * @param  potential The external potential properties.
 */
void potential_print(const struct external_potential* potential) {

#ifdef EXTERNAL_POTENTIAL_POINTMASS
  message("Point mass properties are (x,y,z) = (%e, %e, %e), M = %e",
          potential->point_mass.x, potential->point_mass.y,
          potential->point_mass.z, potential->point_mass.mass);
#endif /* EXTERNAL_POTENTIAL_POINTMASS */
}