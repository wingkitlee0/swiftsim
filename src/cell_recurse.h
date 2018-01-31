/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *               2016 John A. Regan (john.a.regan@durham.ac.uk)
 *                    Tom Theuns (tom.theuns@durham.ac.uk)
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
#ifndef SWIFT_CELL_RECURSE_H
#define SWIFT_CELL_RECURSE_H

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "active.h"
#include "cell.h"
#include "engine.h"

/* Function type for mappings. */
typedef void (*cell_pair_map_function)(struct cell *ci, struct cell *cj,
                                       int sid, void *data);

/**
 * @brief Map a function to all the active cell pairs recursively.
 *
 * @param ci The first #cell in the pair.
 * @param cj The second #cell in the pair.
 * @param e The engine in which these cells exist.
 * @param fun The function to map onto active pairs.
 * @param data Extra data to pass to the mapping function.
 */
static void cell_active_hydro_pairs_recurse(struct cell *ci, struct cell *cj,
                                            const struct engine *e,
                                            cell_pair_map_function fun,
                                            void *data) {

  /* Should we even bother? */
  if (!cell_is_active_hydro(ci, e) && !cell_is_active_hydro(cj, e)) return;
  if (ci->count == 0 || cj->count == 0) return;

  /* Get the type of pair if not specified explicitly. */
  double shift[3];
  const int sid = space_getsid(e->s, &ci, &cj, shift);

  /* Recurse? */
  if (cell_can_recurse_in_pair_task(ci) && cell_can_recurse_in_pair_task(cj)) {

    /* Different types of flags. */
    switch (sid) {

      /* Regular sub-cell interactions of a single cell. */
      case 0: /* (  1 ,  1 ,  1 ) */
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[0], e,
                                          fun, data);
        break;

      case 1: /* (  1 ,  1 ,  0 ) */
        if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[0], e,
                                          fun, data);
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[1], e,
                                          fun, data);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[0], e,
                                          fun, data);
        if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[1], e,
                                          fun, data);
        break;

      case 2: /* (  1 ,  1 , -1 ) */
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[1], e,
                                          fun, data);
        break;

      case 3: /* (  1 ,  0 ,  1 ) */
        if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[5], cj->progeny[0], e,
                                          fun, data);
        if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[5], cj->progeny[2], e,
                                          fun, data);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[0], e,
                                          fun, data);
        if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[2], e,
                                          fun, data);
        break;

      case 4: /* (  1 ,  0 ,  0 ) */
        if (ci->progeny[4] != NULL && cj->progeny[0] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[4], cj->progeny[0], e,
                                          fun, data);
        if (ci->progeny[4] != NULL && cj->progeny[1] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[4], cj->progeny[1], e,
                                          fun, data);
        if (ci->progeny[4] != NULL && cj->progeny[2] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[4], cj->progeny[2], e,
                                          fun, data);
        if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[4], cj->progeny[3], e,
                                          fun, data);
        if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[5], cj->progeny[0], e,
                                          fun, data);
        if (ci->progeny[5] != NULL && cj->progeny[1] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[5], cj->progeny[1], e,
                                          fun, data);
        if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[5], cj->progeny[2], e,
                                          fun, data);
        if (ci->progeny[5] != NULL && cj->progeny[3] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[5], cj->progeny[3], e,
                                          fun, data);
        if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[0], e,
                                          fun, data);
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[1], e,
                                          fun, data);
        if (ci->progeny[6] != NULL && cj->progeny[2] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[2], e,
                                          fun, data);
        if (ci->progeny[6] != NULL && cj->progeny[3] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[3], e,
                                          fun, data);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[0], e,
                                          fun, data);
        if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[1], e,
                                          fun, data);
        if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[2], e,
                                          fun, data);
        if (ci->progeny[7] != NULL && cj->progeny[3] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[3], e,
                                          fun, data);
        break;

      case 5: /* (  1 ,  0 , -1 ) */
        if (ci->progeny[4] != NULL && cj->progeny[1] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[4], cj->progeny[1], e,
                                          fun, data);
        if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[4], cj->progeny[3], e,
                                          fun, data);
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[1], e,
                                          fun, data);
        if (ci->progeny[6] != NULL && cj->progeny[3] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[3], e,
                                          fun, data);
        break;

      case 6: /* (  1 , -1 ,  1 ) */
        if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[5], cj->progeny[2], e,
                                          fun, data);
        break;

      case 7: /* (  1 , -1 ,  0 ) */
        if (ci->progeny[4] != NULL && cj->progeny[2] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[4], cj->progeny[2], e,
                                          fun, data);
        if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[4], cj->progeny[3], e,
                                          fun, data);
        if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[5], cj->progeny[2], e,
                                          fun, data);
        if (ci->progeny[5] != NULL && cj->progeny[3] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[5], cj->progeny[3], e,
                                          fun, data);
        break;

      case 8: /* (  1 , -1 , -1 ) */
        if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[4], cj->progeny[3], e,
                                          fun, data);
        break;

      case 9: /* (  0 ,  1 ,  1 ) */
        if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[3], cj->progeny[0], e,
                                          fun, data);
        if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[3], cj->progeny[4], e,
                                          fun, data);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[0], e,
                                          fun, data);
        if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[4], e,
                                          fun, data);
        break;

      case 10: /* (  0 ,  1 ,  0 ) */
        if (ci->progeny[2] != NULL && cj->progeny[0] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[2], cj->progeny[0], e,
                                          fun, data);
        if (ci->progeny[2] != NULL && cj->progeny[1] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[2], cj->progeny[1], e,
                                          fun, data);
        if (ci->progeny[2] != NULL && cj->progeny[4] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[2], cj->progeny[4], e,
                                          fun, data);
        if (ci->progeny[2] != NULL && cj->progeny[5] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[2], cj->progeny[5], e,
                                          fun, data);
        if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[3], cj->progeny[0], e,
                                          fun, data);
        if (ci->progeny[3] != NULL && cj->progeny[1] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[3], cj->progeny[1], e,
                                          fun, data);
        if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[3], cj->progeny[4], e,
                                          fun, data);
        if (ci->progeny[3] != NULL && cj->progeny[5] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[3], cj->progeny[5], e,
                                          fun, data);
        if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[0], e,
                                          fun, data);
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[1], e,
                                          fun, data);
        if (ci->progeny[6] != NULL && cj->progeny[4] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[4], e,
                                          fun, data);
        if (ci->progeny[6] != NULL && cj->progeny[5] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[5], e,
                                          fun, data);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[0], e,
                                          fun, data);
        if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[1], e,
                                          fun, data);
        if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[4], e,
                                          fun, data);
        if (ci->progeny[7] != NULL && cj->progeny[5] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[5], e,
                                          fun, data);
        break;

      case 11: /* (  0 ,  1 , -1 ) */
        if (ci->progeny[2] != NULL && cj->progeny[1] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[2], cj->progeny[1], e,
                                          fun, data);
        if (ci->progeny[2] != NULL && cj->progeny[5] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[2], cj->progeny[5], e,
                                          fun, data);
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[1], e,
                                          fun, data);
        if (ci->progeny[6] != NULL && cj->progeny[5] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[5], e,
                                          fun, data);
        break;

      case 12: /* (  0 ,  0 ,  1 ) */
        if (ci->progeny[1] != NULL && cj->progeny[0] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[1], cj->progeny[0], e,
                                          fun, data);
        if (ci->progeny[1] != NULL && cj->progeny[2] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[1], cj->progeny[2], e,
                                          fun, data);
        if (ci->progeny[1] != NULL && cj->progeny[4] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[1], cj->progeny[4], e,
                                          fun, data);
        if (ci->progeny[1] != NULL && cj->progeny[6] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[1], cj->progeny[6], e,
                                          fun, data);
        if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[3], cj->progeny[0], e,
                                          fun, data);
        if (ci->progeny[3] != NULL && cj->progeny[2] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[3], cj->progeny[2], e,
                                          fun, data);
        if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[3], cj->progeny[4], e,
                                          fun, data);
        if (ci->progeny[3] != NULL && cj->progeny[6] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[3], cj->progeny[6], e,
                                          fun, data);
        if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[5], cj->progeny[0], e,
                                          fun, data);
        if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[5], cj->progeny[2], e,
                                          fun, data);
        if (ci->progeny[5] != NULL && cj->progeny[4] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[5], cj->progeny[4], e,
                                          fun, data);
        if (ci->progeny[5] != NULL && cj->progeny[6] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[5], cj->progeny[6], e,
                                          fun, data);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[0], e,
                                          fun, data);
        if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[2], e,
                                          fun, data);
        if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[4], e,
                                          fun, data);
        if (ci->progeny[7] != NULL && cj->progeny[6] != NULL)
          cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[6], e,
                                          fun, data);
        break;
    }
  }

  /* Otherwise, compute the pair directly. */
  else if (cell_is_active_hydro(ci, e) || cell_is_active_hydro(cj, e)) {

    /* Call the function on the cell pair. */
    fun(ci, cj, sid, data);
  }
}

#endif /* SWIFT_CELL_RECURSE_H */
