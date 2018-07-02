/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_COOLING_RATES_EAGLE_H
#define SWIFT_COOLING_RATES_EAGLE_H

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "chemistry_struct.h"
#include "cooling_struct.h"

double eagle_do_cooling(const struct cooling_function_data* cooling,
                        const double uold_cgs, const double rho_cgs,
                        const double dt_cgs, const double delta_z,
                        const double redshift,
                        const float Z[chemistry_element_count]);

#endif /*SWIFT_COOLING_RATES_EAGLE_H */
