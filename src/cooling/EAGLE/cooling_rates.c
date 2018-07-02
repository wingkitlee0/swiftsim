
/* Local headers */
#include "cooling.h"
#include "cooling_interpolation.h"
#include "cooling_tables.h"

/*! Relative change in internal energy below which we allow the use of the
 * explicit time-integration solution */
#define eagle_cooling_explicit_tolerance 0.05f

/*! Relative change in internal energy we want the scheme to converge to */
#define eagle_cooling_implicit_tolerance 1.e-6f

/*! Factor used to increase/decrease the bracketing values in the implicit
 * solver */
#define eagle_cooling_bracketing_factor 1.1f

/**
 * @brief Computes the extra heat from Helium reionisation at a given redshift.
 *
 * We follow the implementation of Wiersma et al. 2009, MNRAS, 399, 574-600,
 * section. 2. The calculation returns energy in CGS.
 *
 * Note that delta_z is negative.
 *
 * @param z The current redshift.
 * @param delta_z The change in redhsift over the course of this time-step.
 * @param cooling The #cooling_function_data used in the run.
 */
INLINE static float eagle_cooling_helium_reion_extra_heat(
    const float z, const float delta_z,
    const struct cooling_function_data* restrict cooling) {

#ifdef SWIFT_DEBUG_CHECKS
  if (delta_z > 0.f) error("Invalid value for delta_z. Should be negative.");
#endif

  /* Recover the values we need */
  const float He_reion_z_centre = cooling->He_reion_z_centre;
  const float He_reion_z_sigma = cooling->He_reion_z_sigma;
  const float He_reion_heat_cgs = cooling->He_reion_heat_cgs;

  // MATTHIEU: to do: Optimize this.

  float extra_heat;

  /* Integral of the Gaussian between z and z - delta_z */
  extra_heat =
      erff((z - delta_z - He_reion_z_centre) / (M_SQRT2 * He_reion_z_sigma));
  extra_heat -= erff((z - He_reion_z_centre) / (M_SQRT2 * He_reion_z_sigma));

  /* Multiply by the normalisation factor */
  extra_heat *= He_reion_heat_cgs * 0.5;

  return extra_heat;
}

/**
 * @brief Computes the temperature corresponding to a given internal energy,
 * hydrogen number density, Helium fraction and redshift.
 *
 * Note that the redshift is implicitly passed in via the currently loaded
 * tables in the #cooling_function_data.
 *
 * We interpolate the flattened 4D table 'u_to_temp' that is arranged in the
 * following way:
 * - 1st dim: redshift, length = 2
 * - 2nd dim: Helium fraction, length = eagle_cooling_N_He_frac
 * - 3rd dim: Hydrogen density, length = eagle_cooling_N_density
 * - 4th dim: Internal energy, length = eagle_cooling_N_temperature
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param u_cgs The internal energy in CGS units.
 * @param n_H_cgs The Hydrogen number density in CGS units.
 * @param He_frac The Helium fraction.
 */
float eagle_convert_u_to_T(const struct cooling_function_data* cooling,
                           const float u_cgs, const float n_H_cgs,
                           const float He_frac) {

  /* Start by recovering indices along the axis of the table */
  int index_He, index_nH, index_u;
  float delta_He_table, delta_nH_table, delta_u_table;

  /* Helium fraction axis */
  find_1d_index(cooling->table_He_frac, eagle_cooling_N_He_frac, He_frac,
                &index_He, &delta_He_table);

  /* Hydrogen number density axis */
  find_1d_index(cooling->table_nH, eagle_cooling_N_density, log10f(n_H_cgs),
                &index_nH, &delta_nH_table);

  /* Internal energy axis */
  find_1d_index(cooling->table_u, eagle_cooling_N_temperature, log10f(u_cgs),
                &index_u, &delta_u_table);

  /* For the redshift, we can use the pre-computed delta stored in the cooling
   * function */
  float delta_z_table = cooling->delta_z_table;

  /* Now, interpolate ! */

  /* Note: The 2 for redshifts is because we only have 2 tables in memory at a
   * given point in time */
  const float logT =
      interpolation_4d(cooling->table_u_to_temp, 0, index_He, index_nH, index_u,
                       2, eagle_cooling_N_He_frac, eagle_cooling_N_density,
                       eagle_cooling_N_temperature, delta_z_table,
                       delta_He_table, delta_nH_table, delta_u_table);

  const float T = powf(10.f, logT);

  /* Special case for the low-energy range */
  if (index_u == 0 && delta_u_table == 0.f)
    return T * u_cgs / powf(10.f, cooling->table_u[0]);
  else
    return T;
}

/**
 * @brief Computes the cooling corresponding to a given internal energy,
 * hydrogen number density, Helium fraction, redshift and metallicity.
 *
 * We interpolate the flattened 4D table 'u_to_temp' that is arranged in the
 * following way:
 * - 1st dim: redshift, length = 2
 * - 2nd dim: Helium fraction, length = eagle_cooling_N_He_frac
 * - 3rd dim: Hydrogen density, length = eagle_cooling_N_density
 * - 4th dim: Internal energy, length = eagle_cooling_N_temperature
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param u_cgs The internal energy in CGS units.
 * @param n_H_cgs The Hydrogen number density in CGS units.
 * @param He_frac The Helium fraction.
 * @param redshift The current redshift.
 * @param delta_z The time-step expressed as a change in redshift. This is
 * negative.
 * @param Z The element metallicity expressed as element mass fractions.
 */
float eagle_cooling_rate(const struct cooling_function_data* cooling,
                         const float u_cgs, const float n_H_cgs,
                         const float He_frac, const float redshift,
                         const float delta_z,
                         const float Z[chemistry_element_count]) {

  /* Get the temperature corresponding to this energy */
  const float temp = eagle_convert_u_to_T(cooling, u_cgs, n_H_cgs, He_frac);

  /* Start by recovering indices along the axis of the tables */
  int index_T, index_He, index_nH;
  float delta_T_table, delta_He_table, delta_nH_table;

  /* Temperature axis */
  find_1d_index(cooling->table_temperatures, eagle_cooling_N_temperature,
                log10f(temp), &index_T, &delta_T_table);

  /* Helium fraction axis */
  find_1d_index(cooling->table_He_frac, eagle_cooling_N_He_frac, He_frac,
                &index_He, &delta_He_table);

  /* Hydrogen number density axis */
  find_1d_index(cooling->table_nH, eagle_cooling_N_density, log10f(n_H_cgs),
                &index_nH, &delta_nH_table);

  /* For the redshift, we can use the pre-computed delta stored in the cooling
   * function */
  float delta_z_table = cooling->delta_z_table;

  /**********************/
  /* Metal-free cooling */
  /**********************/

  float Lambda_net = interpolation_4d(
      cooling->table_H_and_He_net_heating, 0, index_He, index_nH, index_T, 2,
      eagle_cooling_N_He_frac, eagle_cooling_N_density,
      eagle_cooling_N_temperature, delta_z_table, delta_He_table,
      delta_nH_table, delta_T_table);

  /**********************/
  /* Compton cooling    */
  /**********************/

  /* Do we need to add the inverse Compton cooling? */
  /* It is *not* stored in the tables before re-ionisation */
  if ((redshift > cooling->table_redshifts[eagle_cooling_N_redshifts - 1]) ||
      (redshift > cooling->H_reion_z)) {

    const float electron_abundance = interpolation_4d(
        cooling->table_H_and_He_electron_abundance, 0, index_He, index_nH,
        index_T, 2, eagle_cooling_N_He_frac, eagle_cooling_N_density,
        eagle_cooling_N_temperature, delta_z_table, delta_He_table,
        delta_nH_table, delta_T_table);

    const float zp1 = 1.f + redshift;
    const float zp1p4 = zp1 * zp1 * zp1 * zp1;

    /* CMB temperature at this redshift */
    const float T_CMB = cooling->T_CMB_0 * zp1;

    /* Compton cooling rate */
    Lambda_net -= cooling->compton_rate_cgs * (temp - T_CMB) * zp1p4 *
                  electron_abundance / n_H_cgs;
  }

  /**********************/
  /* Metal-line cooling */
  /**********************/

  return Lambda_net;
}

/**
 * @brief Compute the new energy of a particle based on its current
 * properties.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param uold_cgs The old internal energy per unit mass in CGS units.
 * @param rho_cgs The current density in CGS units.
 * @param dt_cgs The time-step over which to cool in CGS units.
 * @param delta_z The time-step expressed as a change in redshift. This is
 * negative.
 * @param redshift The current redshift.
 * @param Z The element metallicity expressed as element mass fractions.
 *
 * @return The new internal energy in CGS units.
 */
float eagle_do_cooling(const struct cooling_function_data* cooling,
                       const float uold_cgs, const float rho_cgs,
                       const float dt_cgs, const float delta_z,
                       const float redshift,
                       const float Z[chemistry_element_count]) {

#ifdef SWIFT_DEBUG_CHECKS
  if (delta_z > 0.f) error("Invalid value for delta_z. Should be negative.");
#endif

  /* Hydrogen fraction */
  const float XH = Z[chemistry_element_H];

  /* Helium fraction */
  const float He_frac =
      Z[chemistry_element_He] / (XH + Z[chemistry_element_He]);

  /* Hydrogen number density */
  const float n_H_cgs = rho_cgs * XH / cooling->const_proton_mass_cgs;

  /* Cooling rate factor n_H^2 / rho = n_H * XH / m_p */
  const float rate_factor_cgs = n_H_cgs * XH / cooling->const_proton_mass_cgs;

  /* Compute the extra heat injected by Helium re-ionization */
  const float u_reion_cgs =
      eagle_cooling_helium_reion_extra_heat(redshift, delta_z, cooling);

  /* Turn this into a cooling (heating) rate */
  const float Lambda_reion = u_reion_cgs / (dt_cgs * rate_factor_cgs);

  /*************************************/
  /* Let's get a first guess           */
  /*************************************/

  /* First guess */
  float Lambda_net =
      Lambda_reion + eagle_cooling_rate(cooling, uold_cgs, n_H_cgs, He_frac,
                                        redshift, delta_z, Z);

  float delta_u_cgs = rate_factor_cgs * Lambda_net * dt_cgs;

  /* Small cooling rate: Use explicit solution and abort */
  if (fabsf(delta_u_cgs) < eagle_cooling_explicit_tolerance * uold_cgs) {

    /* Explicit solution */
    return uold_cgs + delta_u_cgs;
  }

  /*************************************/
  /* Let's try to bracket the solution */
  /*************************************/

  float u_cgs = uold_cgs;
  float u_lower_cgs = u_cgs;
  float u_upper_cgs = u_cgs;

  if (delta_u_cgs > 0.f) { /* Net heating case */

    u_upper_cgs *= sqrtf(eagle_cooling_bracketing_factor);
    u_lower_cgs /= sqrtf(eagle_cooling_bracketing_factor);

    /* Compute a new rate */
    float Lambda_new = eagle_cooling_rate(cooling, u_upper_cgs, n_H_cgs,
                                          He_frac, redshift, delta_z, Z);
    float delta_new_cgs = u_reion_cgs + Lambda_new * rate_factor_cgs * dt_cgs;

    while (u_upper_cgs - uold_cgs - delta_new_cgs < 0.f) {

      u_upper_cgs *= eagle_cooling_bracketing_factor;
      u_lower_cgs *= eagle_cooling_bracketing_factor;

      /* Compute a new rate */
      Lambda_new = eagle_cooling_rate(cooling, u_upper_cgs, n_H_cgs, He_frac,
                                      redshift, delta_z, Z);
      delta_new_cgs = u_reion_cgs + Lambda_new * rate_factor_cgs * dt_cgs;
    }
  }

  if (delta_u_cgs < 0.f) { /* Net cooling case */

    u_upper_cgs *= sqrtf(eagle_cooling_bracketing_factor);
    u_lower_cgs /= sqrtf(eagle_cooling_bracketing_factor);

    /* Compute a new rate */
    float Lambda_new = eagle_cooling_rate(cooling, u_lower_cgs, n_H_cgs,
                                          He_frac, redshift, delta_z, Z);
    float delta_new_cgs = u_reion_cgs + Lambda_new * rate_factor_cgs * dt_cgs;

    while (u_lower_cgs - uold_cgs - delta_new_cgs > 0.f) {

      u_upper_cgs /= eagle_cooling_bracketing_factor;
      u_lower_cgs /= eagle_cooling_bracketing_factor;

      /* Compute a new rate */
      Lambda_new = eagle_cooling_rate(cooling, u_lower_cgs, n_H_cgs, He_frac,
                                      redshift, delta_z, Z);
      delta_new_cgs = u_reion_cgs + Lambda_new * rate_factor_cgs * dt_cgs;
    }
  }

  /********************************************/
  /* We now have an upper and lower bound.    */
  /* Let's iterate by reducing the bracketing */
  /********************************************/

  float delta_upper_lower = 0.f;
  do {

    /* New guess */
    u_cgs = 0.5f * (u_lower_cgs + u_upper_cgs);

    /* New rate */
    Lambda_net =
        Lambda_reion + eagle_cooling_rate(cooling, u_cgs, n_H_cgs, He_frac,
                                          redshift, delta_z, Z);

    /* New change in energy */
    delta_u_cgs = rate_factor_cgs * Lambda_net * dt_cgs;

    /* Bracketing */
    if (u_cgs - uold_cgs - delta_u_cgs > 0.f)
      u_upper_cgs = u_cgs;
    else
      u_lower_cgs = u_cgs;

    /* Width of the range */
    delta_upper_lower = u_upper_cgs - u_lower_cgs;

  } while (fabsf(delta_upper_lower) > eagle_cooling_implicit_tolerance * u_cgs);

  return u_cgs;
}
