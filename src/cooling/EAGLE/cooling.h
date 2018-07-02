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
#ifndef SWIFT_COOLING_EAGLE_H
#define SWIFT_COOLING_EAGLE_H

/**
 * @file src/cooling/EAGLE/cooling.h
 * @brief EAGLE cooling function
 */

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <float.h>
#include <math.h>

/* Local includes. */
#include "cooling_rates.h"
#include "cooling_tables.h"
#include "error.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

/**
 * @brief Common operations performed on the cooling function at a
 * given time-step or redshift.
 *
 * Here we load the cooling tables corresponding to our redshift.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 */
INLINE static void cooling_update(const struct phys_const* phys_const,
                                  const struct unit_system* us,
                                  const struct cosmology* cosmo,
                                  struct cooling_function_data* cooling) {

  /* Current redshift */
  const float redshift = cosmo->z;

  /* Get index along the redshift index of the table */
  int index_z;
  float delta_z_table;
  eagle_get_redshift_table_index(redshift, cooling, &index_z, &delta_z_table);
  cooling->index_z = index_z;
  cooling->delta_z_table = delta_z_table;

  /* Load the current table (if different from what we have) */
  ealge_check_cooling_tables(cooling, index_z);
}

/**
 * @brief Apply the cooling function to a particle.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param dt The time-step of this particle (in internal units).
 */
__attribute__((always_inline)) INLINE static void cooling_cool_part(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct cooling_function_data* restrict cooling,
    struct part* restrict p, struct xpart* restrict xp, float dt) {

  /* Abort early if the time-step is 0 */
  if (dt == 0.f) return;

  /* Current redshift and change in redshift of this time-step */
  const float redshift = cosmo->z;
  const float delta_z = -dt / cosmo->a;  // MATTHIEU: Check this!

  /* Recover some particle properties */
  const float mass = hydro_get_mass(p);
  const float mass_inv = 1.f / mass;
  const float rho = hydro_get_physical_density(p, cosmo);
  const float uold = hydro_get_physical_internal_energy(p, cosmo);

  /* Conversion to CGS system */
  const float dt_cgs = dt * units_cgs_conversion_factor(us, UNIT_CONV_TIME);
  const float rho_cgs =
      rho * units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
  const float uold_cgs =
      uold * units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);

  // MATTHIEU: to do: Add check for min energy here

  /* Compute metallicites as mass fractions of the smoothed metallicites */
  float Z[chemistry_element_count];
  for (int i = 0; i < chemistry_element_count; ++i)
    Z[i] = p->chemistry_data.smoothed_metal_mass_fraction[chemistry_element_H] *
           mass_inv;

  /* Do the cooling */
  const float unew_cgs = eagle_do_cooling(cooling, uold_cgs, rho_cgs, dt_cgs,
                                          delta_z, redshift, Z);

  // MATTHIEU: to do: Add check for min energy here

  /* Write back to the particle */
  const float unew = unew_cgs / units_cgs_conversion_factor(
                                    us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  hydro_set_internal_energy_dt(p, (unew - uold) / dt);
}

/**
 * @brief Computes the cooling time-step.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float cooling_timestep(
    const struct cooling_function_data* restrict cooling,
    const struct phys_const* restrict phys_const,
    const struct cosmology* restrict cosmo,
    const struct unit_system* restrict us, const struct part* restrict p) {

  return FLT_MAX;
}

/**
 * @brief Sets the cooling properties of the (x-)particles to a valid start
 * state.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void cooling_first_init_part(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct cooling_function_data* restrict cooling,
    const struct part* restrict p, struct xpart* restrict xp) {}

/**
 * @brief Returns the total radiated energy by this particle.
 *
 * @param xp The extended particle data
 */
__attribute__((always_inline)) INLINE static float cooling_get_radiated_energy(
    const struct xpart* restrict xp) {

  return 0.f;
}

/**
 * @brief Initialises the cooling properties.
 *
 * @param parameter_file The parsed parameter file.
 * @param us The current internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param cooling The cooling properties to initialize
 */
static INLINE void cooling_init_backend(struct swift_params* parameter_file,
                                        const struct unit_system* us,
                                        const struct phys_const* phys_const,
                                        struct cooling_function_data* cooling) {

  /* Get the proton mass in CGS */
  cooling->const_proton_mass_cgs =
      phys_const->const_proton_mass *
      units_cgs_conversion_factor(us, UNIT_CONV_MASS);

  /* Read parameters related to H reionisation */
  cooling->H_reion_z =
      parser_get_param_double(parameter_file, "EagleCooling:H_reion_z");

  /* Read parameters related to He II reionisation */
  cooling->He_reion_z_centre =
      parser_get_param_double(parameter_file, "EagleCooling:He_reion_z_centre");
  cooling->He_reion_z_sigma =
      parser_get_param_double(parameter_file, "EagleCooling:He_reion_z_sigma");
  cooling->He_reion_heat_cgs = parser_get_param_double(
      parameter_file, "EagleCooling:He_reion_heat_eVpH");

  /* Convert He II reionisation parameters to the units we use internally */
  cooling->He_reion_heat_cgs *=
      phys_const->const_electron_volt / phys_const->const_proton_mass;
  cooling->He_reion_heat_cgs *=
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);

  /* Read the cooling tables */
  eagle_cooling_init_redshift_tables(cooling);
  eagle_read_cooling_table_header(cooling);
  eagle_set_solar_metallicity(cooling);
  eagle_allocate_cooling_tables(cooling);

  /* Set the table parameters to un-initialised values */
  cooling->index_z = -1;
  cooling->delta_z_table = 0.f;
  cooling->low_z_index = -1;
  cooling->high_z_index = -1;
}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * @param cooling The properties of the cooling function.
 */
static INLINE void cooling_print_backend(
    const struct cooling_function_data* cooling) {

  message("Cooling function is 'EAGLE'.");

  message("Helium reionisation centre: z=%.2f sigma=%.2f",
          cooling->He_reion_z_centre, cooling->He_reion_z_sigma);
  message("Helium reionisation heating: %e [erg/g]",
          cooling->He_reion_heat_cgs);
}

#endif /* SWIFT_COOLING_EAGLE_H */
