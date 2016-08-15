/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016   Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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
#ifndef SWIFT_EQUATION_OF_STATE_H
#define SWIFT_EQUATION_OF_STATE_H

/**
 * @file equation_of_state.h
 * @brief Defines the equation of state of the gas we simulate in the form of
 * relations between thermodynamic quantities. These are later used internally
 * by all hydro schemes
 */

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <math.h>

/* Local headers. */
#include "adiabatic_index.h"
#include "const.h"
#include "debug.h"
#include "inline.h"

/* ------------------------------------------------------------------------- */
#if defined(EOS_IDEAL_GAS)

/**
 * @brief Returns the internal energy given density and entropy
 *
 * Computes \f$u = \frac{S\rho^{\gamma-1} }{\gamma - 1}\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$S\f$.
 */
__attribute__((always_inline)) INLINE static float
gas_internal_energy_from_entropy(float density, float entropy) {

  return entropy * pow_gamma_minus_one(density) *
         hydro_one_over_gamma_minus_one;
}

/**
 * @brief Returns the pressure given density and entropy
 *
 * Computes \f$P = S\rho^\gamma\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$S\f$.
 */
__attribute__((always_inline)) INLINE static float gas_pressure_from_entropy(
    float density, float entropy) {

  return entropy * pow_gamma(density);
}

/**
 * @brief Returns the sound speed given density and entropy
 *
 * Computes \f$c = \sqrt{\gamma S \rho^{\gamma-1}}\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$S\f$.
 */
__attribute__((always_inline)) INLINE static float gas_soundspeed_from_entropy(
    float density, float entropy) {

  return sqrtf(hydro_gamma * pow_gamma_minus_one(density) * entropy);
}

/**
 * @brief Returns the entropy given density and internal energy
 *
 * Computes \f$S = \frac{(\gamma - 1)u}{\rho^{\gamma-1}}\f$.
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline)) INLINE static float
gas_entropy_from_internal_energy(float density, float u) {

  return hydro_gamma_minus_one * u * pow_minus_gamma_minus_one(density);
}

/**
 * @brief Returns the pressure given density and internal energy
 *
 * Computes \f$P = (\gamma - 1)u\rho\f$.
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline)) INLINE static float
gas_pressure_from_internal_energy(float density, float u) {

  return hydro_gamma_minus_one * u * density;
}

/**
 * @brief Returns the sound speed given density and internal energy
 *
 * Computes \f$c = \sqrt{\gamma (\gamma - 1) u }\f$.
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline)) INLINE static float
gas_soundspeed_from_internal_energy(float density, float u) {

  return sqrtf(u * hydro_gamma * hydro_gamma_minus_one);
}

/* ------------------------------------------------------------------------- */
#elif defined(EOS_ISOTHERMAL_GAS)

/**
 * @brief Returns the internal energy given density and entropy
 *
 * @param density The density
 * @param entropy The entropy
 */
__attribute__((always_inline)) INLINE static float
gas_internal_energy_from_entropy(float density, float entropy) {

  error("Missing definition !");
  return 0.f;
}

/**
 * @brief Returns the pressure given density and entropy
 *
 * @param density The density
 * @param entropy The entropy
 */
__attribute__((always_inline)) INLINE static float gas_pressure_from_entropy(
    float density, float entropy) {

  error("Missing definition !");
  return 0.f;
}

/**
 * @brief Returns the entropy given density and internal energy
 *
 * @param density The density
 * @param u The internal energy
 */
__attribute__((always_inline)) INLINE static float
gas_entropy_from_internal_energy(float density, float u) {

  error("Missing definition !");
  return 0.f;
}

/**
 * @brief Returns the pressure given density and internal energy
 *
 * @param density The density
 * @param u The internal energy
 */
__attribute__((always_inline)) INLINE static float
gas_pressure_from_internal_energy(float density, float u) {

  error("Missing definition !");
  return 0.f;
}

/**
 * @brief Returns the sound speed given density and internal energy
 *
 * @param density The density
 * @param u The internal energy
 */
__attribute__((always_inline)) INLINE static float
gas_soundspeed_from_internal_energy(float density, float u) {

  error("Missing definition !");
  return 0.f;
}

/* ------------------------------------------------------------------------- */
#else

#error "An Equation of state needs to be chosen in const.h !"

#endif

#endif /* SWIFT_EQUATION_OF_STATE_H */