/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_COOLING_TABLES_EAGLE_H
#define SWIFT_COOLING_TABLES_EAGLE_H

#include "cooling_struct.h"

void get_redshift_table_index(
    const float z, const struct cooling_function_data *restrict cooling,
    int *z_index, float *delta_z);

void eagle_cooling_init_redshift_tables(struct cooling_function_data *cooling);

void eagle_read_cooling_table_header(struct cooling_function_data *cooling);

void eagle_set_solar_metallicity(struct cooling_function_data *cooling);

void eagle_allocate_cooling_tables(struct cooling_function_data *cooling);

void ealge_check_cooling_tables(struct cooling_function_data *cooling,
                                int index_z);

#endif /* SWIFT_COOLING_TABLES_EAGLE_H */
