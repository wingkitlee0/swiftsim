/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
 * @param do_self If true, also apply @c fun to each individual cell
 *recursively.
 * @param fun The function to map onto active pairs. If @c do_self is true, may
 *        be called with its second parameter set to @c NULL.
 * @param data Extra data to pass to the mapping function.
 */
static void cell_active_hydro_pairs_recurse(struct cell *ci, struct cell *cj,
                                            const struct engine *e, int do_self,
                                            cell_pair_map_function fun,
                                            void *data) {

  /* Should we even bother? */
  if (!cell_is_active_hydro(ci, e) &&
      (cj == NULL || !cell_is_active_hydro(cj, e)))
    return;
  if (ci->count == 0 || (cj != NULL && cj->count == 0)) return;

  /* Self interaction? */
  if (cj == NULL) {
    /* Do anything? */
    if (!cell_is_active_hydro(ci, e)) return;

    /* Recurse? */
    if (cell_can_recurse_in_self_task(ci)) {

      /* Loop over all progeny and pairs of progenies */
      for (int j = 0; j < 8; j++) {
        if (ci->progeny[j] != NULL) {
          cell_active_hydro_pairs_recurse(ci->progeny[j], NULL, e, do_self, fun,
                                          data);
          for (int k = j + 1; k < 8; k++)
            if (ci->progeny[k] != NULL)
              cell_active_hydro_pairs_recurse(ci->progeny[j], ci->progeny[k], e,
                                              do_self, fun, data);
        }
      }
    } else {

      /* We have reached the bottom of the tree, apply the function. */
      fun(ci, NULL, /* sid */ -1, data);
    }
  } else {

    /* Get the type of pair. */
    double shift[3];
    const int sid = space_getsid(e->s, &ci, &cj, shift);

    /* Recurse? */
    if (cell_can_recurse_in_pair_task(ci) &&
        cell_can_recurse_in_pair_task(cj)) {

      /* Different types of flags. */
      switch (sid) {

        /* Regular sub-cell interactions of a single cell. */
        case 0: /* (  1 ,  1 ,  1 ) */
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[0], e,
                                            do_self, fun, data);
          break;

        case 1: /* (  1 ,  1 ,  0 ) */
          if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[0], e,
                                            do_self, fun, data);
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[1], e,
                                            do_self, fun, data);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[0], e,
                                            do_self, fun, data);
          if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[1], e,
                                            do_self, fun, data);
          break;

        case 2: /* (  1 ,  1 , -1 ) */
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[1], e,
                                            do_self, fun, data);
          break;

        case 3: /* (  1 ,  0 ,  1 ) */
          if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[5], cj->progeny[0], e,
                                            do_self, fun, data);
          if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[5], cj->progeny[2], e,
                                            do_self, fun, data);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[0], e,
                                            do_self, fun, data);
          if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[2], e,
                                            do_self, fun, data);
          break;

        case 4: /* (  1 ,  0 ,  0 ) */
          if (ci->progeny[4] != NULL && cj->progeny[0] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[4], cj->progeny[0], e,
                                            do_self, fun, data);
          if (ci->progeny[4] != NULL && cj->progeny[1] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[4], cj->progeny[1], e,
                                            do_self, fun, data);
          if (ci->progeny[4] != NULL && cj->progeny[2] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[4], cj->progeny[2], e,
                                            do_self, fun, data);
          if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[4], cj->progeny[3], e,
                                            do_self, fun, data);
          if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[5], cj->progeny[0], e,
                                            do_self, fun, data);
          if (ci->progeny[5] != NULL && cj->progeny[1] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[5], cj->progeny[1], e,
                                            do_self, fun, data);
          if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[5], cj->progeny[2], e,
                                            do_self, fun, data);
          if (ci->progeny[5] != NULL && cj->progeny[3] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[5], cj->progeny[3], e,
                                            do_self, fun, data);
          if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[0], e,
                                            do_self, fun, data);
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[1], e,
                                            do_self, fun, data);
          if (ci->progeny[6] != NULL && cj->progeny[2] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[2], e,
                                            do_self, fun, data);
          if (ci->progeny[6] != NULL && cj->progeny[3] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[3], e,
                                            do_self, fun, data);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[0], e,
                                            do_self, fun, data);
          if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[1], e,
                                            do_self, fun, data);
          if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[2], e,
                                            do_self, fun, data);
          if (ci->progeny[7] != NULL && cj->progeny[3] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[3], e,
                                            do_self, fun, data);
          break;

        case 5: /* (  1 ,  0 , -1 ) */
          if (ci->progeny[4] != NULL && cj->progeny[1] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[4], cj->progeny[1], e,
                                            do_self, fun, data);
          if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[4], cj->progeny[3], e,
                                            do_self, fun, data);
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[1], e,
                                            do_self, fun, data);
          if (ci->progeny[6] != NULL && cj->progeny[3] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[3], e,
                                            do_self, fun, data);
          break;

        case 6: /* (  1 , -1 ,  1 ) */
          if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[5], cj->progeny[2], e,
                                            do_self, fun, data);
          break;

        case 7: /* (  1 , -1 ,  0 ) */
          if (ci->progeny[4] != NULL && cj->progeny[2] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[4], cj->progeny[2], e,
                                            do_self, fun, data);
          if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[4], cj->progeny[3], e,
                                            do_self, fun, data);
          if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[5], cj->progeny[2], e,
                                            do_self, fun, data);
          if (ci->progeny[5] != NULL && cj->progeny[3] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[5], cj->progeny[3], e,
                                            do_self, fun, data);
          break;

        case 8: /* (  1 , -1 , -1 ) */
          if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[4], cj->progeny[3], e,
                                            do_self, fun, data);
          break;

        case 9: /* (  0 ,  1 ,  1 ) */
          if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[3], cj->progeny[0], e,
                                            do_self, fun, data);
          if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[3], cj->progeny[4], e,
                                            do_self, fun, data);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[0], e,
                                            do_self, fun, data);
          if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[4], e,
                                            do_self, fun, data);
          break;

        case 10: /* (  0 ,  1 ,  0 ) */
          if (ci->progeny[2] != NULL && cj->progeny[0] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[2], cj->progeny[0], e,
                                            do_self, fun, data);
          if (ci->progeny[2] != NULL && cj->progeny[1] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[2], cj->progeny[1], e,
                                            do_self, fun, data);
          if (ci->progeny[2] != NULL && cj->progeny[4] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[2], cj->progeny[4], e,
                                            do_self, fun, data);
          if (ci->progeny[2] != NULL && cj->progeny[5] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[2], cj->progeny[5], e,
                                            do_self, fun, data);
          if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[3], cj->progeny[0], e,
                                            do_self, fun, data);
          if (ci->progeny[3] != NULL && cj->progeny[1] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[3], cj->progeny[1], e,
                                            do_self, fun, data);
          if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[3], cj->progeny[4], e,
                                            do_self, fun, data);
          if (ci->progeny[3] != NULL && cj->progeny[5] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[3], cj->progeny[5], e,
                                            do_self, fun, data);
          if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[0], e,
                                            do_self, fun, data);
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[1], e,
                                            do_self, fun, data);
          if (ci->progeny[6] != NULL && cj->progeny[4] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[4], e,
                                            do_self, fun, data);
          if (ci->progeny[6] != NULL && cj->progeny[5] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[5], e,
                                            do_self, fun, data);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[0], e,
                                            do_self, fun, data);
          if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[1], e,
                                            do_self, fun, data);
          if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[4], e,
                                            do_self, fun, data);
          if (ci->progeny[7] != NULL && cj->progeny[5] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[5], e,
                                            do_self, fun, data);
          break;

        case 11: /* (  0 ,  1 , -1 ) */
          if (ci->progeny[2] != NULL && cj->progeny[1] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[2], cj->progeny[1], e,
                                            do_self, fun, data);
          if (ci->progeny[2] != NULL && cj->progeny[5] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[2], cj->progeny[5], e,
                                            do_self, fun, data);
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[1], e,
                                            do_self, fun, data);
          if (ci->progeny[6] != NULL && cj->progeny[5] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[6], cj->progeny[5], e,
                                            do_self, fun, data);
          break;

        case 12: /* (  0 ,  0 ,  1 ) */
          if (ci->progeny[1] != NULL && cj->progeny[0] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[1], cj->progeny[0], e,
                                            do_self, fun, data);
          if (ci->progeny[1] != NULL && cj->progeny[2] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[1], cj->progeny[2], e,
                                            do_self, fun, data);
          if (ci->progeny[1] != NULL && cj->progeny[4] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[1], cj->progeny[4], e,
                                            do_self, fun, data);
          if (ci->progeny[1] != NULL && cj->progeny[6] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[1], cj->progeny[6], e,
                                            do_self, fun, data);
          if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[3], cj->progeny[0], e,
                                            do_self, fun, data);
          if (ci->progeny[3] != NULL && cj->progeny[2] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[3], cj->progeny[2], e,
                                            do_self, fun, data);
          if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[3], cj->progeny[4], e,
                                            do_self, fun, data);
          if (ci->progeny[3] != NULL && cj->progeny[6] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[3], cj->progeny[6], e,
                                            do_self, fun, data);
          if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[5], cj->progeny[0], e,
                                            do_self, fun, data);
          if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[5], cj->progeny[2], e,
                                            do_self, fun, data);
          if (ci->progeny[5] != NULL && cj->progeny[4] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[5], cj->progeny[4], e,
                                            do_self, fun, data);
          if (ci->progeny[5] != NULL && cj->progeny[6] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[5], cj->progeny[6], e,
                                            do_self, fun, data);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[0], e,
                                            do_self, fun, data);
          if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[2], e,
                                            do_self, fun, data);
          if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[4], e,
                                            do_self, fun, data);
          if (ci->progeny[7] != NULL && cj->progeny[6] != NULL)
            cell_active_hydro_pairs_recurse(ci->progeny[7], cj->progeny[6], e,
                                            do_self, fun, data);
          break;
      }
    }

    /* Otherwise, compute the pair directly. */
    else if (cell_is_active_hydro(ci, e) || cell_is_active_hydro(cj, e)) {

      /* Call the function on the cell pair. */
      fun(ci, cj, sid, data);
    }
  }
}

#endif /* SWIFT_CELL_RECURSE_H */
