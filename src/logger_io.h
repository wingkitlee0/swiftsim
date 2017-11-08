/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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
#ifndef SWIFT_LOGGER_IO_H
#define SWIFT_LOGGER_IO_H

/* Config parameters. */
#include "../config.h"

#if defined(HAVE_HDF5) && !defined(WITH_MPI) && defined(WITH_LOGGER)

/* Includes. */
#include "engine.h"
#include "part.h"
#include "units.h"

void read_ic_single(char* fileName, const struct unit_system* internal_units,
                    double dim[3], struct part** parts, struct gpart** gparts,
                    struct spart** sparts, size_t* Ngas, size_t* Ndm,
                    size_t* Nstars, int* periodic, int* flag_entropy,
                    int with_hydro, int with_gravity, int with_stars,
                    int dry_run);

void write_output_single(struct engine* e, const char* baseName,
                         const struct unit_system* internal_units,
                         const struct unit_system* snapshot_units);

void write_index_single(struct engine* e, const char* baseName,
			const struct unit_system* internal_units,
			const struct unit_system* snapshot_units);
#endif

#endif /* SWIFT_LOGGER_IO_H */