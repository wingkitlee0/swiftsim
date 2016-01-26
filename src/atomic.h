/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#ifndef SWIFT_ATOMIC_H
#define SWIFT_ATOMIC_H

/* Includes. */
#include "inline.h"

#define atomic_add(v, i) __sync_fetch_and_add(v, i)
#define atomic_inc(v) atomic_add(v, 1)
#define atomic_dec(v) atomic_add(v, -1)
#define atomic_cas(v, o, n) __sync_val_compare_and_swap(v, o, n)
#define atomic_load(v) atomic_add(v, 0)

#endif /* SWIFT_ATOMIC_H */
