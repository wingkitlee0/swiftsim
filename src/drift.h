/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_DRIFT_H
#define SWIFT_DRIFT_H

/* Config parameters. */
#include "../config.h"

/* Local headers. */
#include "const.h"
#include "debug.h"
#include "hydro.h"

/**
 * @brief Perform the 'drift' operation on a #gpart
 *
 * @param gp The #gpart to drift.
 * @param dt The drift time-step
 * @param timeBase The minimal allowed time-step size.
 * @param ti_old Integer start of time-step
 * @param ti_current Integer end of time-step
 */
__attribute__((always_inline)) INLINE static void drift_gpart(
    struct gpart* gp, float dt, double timeBase, int ti_old, int ti_current) {
  /* Drift... */
  gp->x[0] += gp->v_full[0] * dt;
  gp->x[1] += gp->v_full[1] * dt;
  gp->x[2] += gp->v_full[2] * dt;

  /* Compute offset since last cell construction */
  gp->x_diff[0] -= gp->v_full[0] * dt;
  gp->x_diff[1] -= gp->v_full[1] * dt;
  gp->x_diff[2] -= gp->v_full[2] * dt;
}

/**
 * @brief Perform the 'drift' operation on a #part
 *
 * @param p The #part to drift.
 * @param xp The #xpart of the particle.
 * @param dt The drift time-step
 * @param timeBase The minimal allowed time-step size.
 * @param ti_old Integer start of time-step
 * @param ti_current Integer end of time-step
 */
__attribute__((always_inline)) INLINE static void drift_part(
    struct part* p, struct xpart* xp, float dt, double timeBase, int ti_old,
    int ti_current) {
  /* Useful quantity */
  const float h_inv = 1.0f / p->h;

  /* Drift... */
  p->x[0] += xp->v_full[0] * dt;
  p->x[1] += xp->v_full[1] * dt;
  p->x[2] += xp->v_full[2] * dt;

  /* Predict velocities (for hydro terms) */
  p->v[0] += p->a_hydro[0] * dt;
  p->v[1] += p->a_hydro[1] * dt;
  p->v[2] += p->a_hydro[2] * dt;

  /* Predict smoothing length */
  const float w1 = p->h_dt * h_inv * dt;
  if (fabsf(w1) < 0.2f)
    p->h *= approx_expf(w1); /* 4th order expansion of exp(w) */
  else
    p->h *= expf(w1);

  /* Predict density */
  const float w2 = -3.0f * w1;
  if (fabsf(w2) < 0.2f)
    p->rho *= approx_expf(w2); /* 4th order expansion of exp(w) */
  else
    p->rho *= expf(w2);

  /* Predict the values of the extra fields */
  hydro_predict_extra(p, xp, ti_old, ti_current, timeBase);

  /* Compute offset since last cell construction */
  xp->x_diff[0] -= xp->v_full[0] * dt;
  xp->x_diff[1] -= xp->v_full[1] * dt;
  xp->x_diff[2] -= xp->v_full[2] * dt;
}

#endif /* SWIFT_DRIFT_H */