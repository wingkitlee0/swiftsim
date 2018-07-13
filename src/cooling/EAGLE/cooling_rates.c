
/* Local headers */
#include "cooling.h"
#include "cooling_interpolation.h"
#include "cooling_tables.h"

/*! Relative change in internal energy below which we allow the use of the
 * explicit time-integration solution */
#define eagle_cooling_explicit_tolerance 0.05

/*! Relative change in internal energy we want the scheme to converge to */
#define eagle_cooling_implicit_tolerance 1.e-6

/*! Factor used to increase/decrease the bracketing values in the implicit
 * solver */
#define eagle_cooling_bracketing_factor 1.1

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
double eagle_cooling_helium_reion_extra_heat(
    const double z, const double delta_z,
    const struct cooling_function_data* restrict cooling) {

#ifdef SWIFT_DEBUG_CHECKS
  if (delta_z > 0.f) error("Invalid value for delta_z. Should be negative.");
#endif

  /* Recover the values we need */
  const double z_centre = cooling->He_reion_z_centre;
  const double z_sigma = cooling->He_reion_z_sigma;
  const double heat_cgs = cooling->He_reion_heat_cgs;

  // MATTHIEU: to do: Optimize this.

  double extra_heat;

  /* Integral of the Gaussian between z and z - delta_z */
  extra_heat = erf((z - delta_z - z_centre) / (M_SQRT2 * z_sigma));
  extra_heat -= erf((z - z_centre) / (M_SQRT2 * z_sigma));

  /* Multiply by the normalisation factor */
  extra_heat *= heat_cgs * 0.5;

  return extra_heat;
}

/**
 * @brief Sets the Metal abundances in units of the solar abundances
 * assumed by the cooling tables.
 *
 * Note that the resulting table does *not* contain H and He since
 * we only want the metals.
 *
 * Calcium and Sulphur are not tracked by the hydro scheme. We use
 * the Silicon abundance for them.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param Z The particle metal content expressed as metal mass fractions.
 * @param element_abundance_solar The abundances in units of solar abundances.
 */
void eagle_set_cooling_abundances(
    const struct cooling_function_data* cooling,
    const float Z[chemistry_element_count],
    float element_abundance_solar[eagle_cooling_N_metal]) {

  /* We copy over the Metals (i.e. ignoring H and He in the first 2 entries
   * of Z[]). For S and Ca we use the Si abundance as we don't track
   * them individually. */

  element_abundance_solar[0] =
      Z[2] / cooling->table_solar_abundances[2]; /* C  */
  element_abundance_solar[1] =
      Z[3] / cooling->table_solar_abundances[3]; /* N  */
  element_abundance_solar[2] =
      Z[4] / cooling->table_solar_abundances[4]; /* O  */
  element_abundance_solar[3] =
      Z[5] / cooling->table_solar_abundances[5]; /* Ne */
  element_abundance_solar[4] =
      Z[6] / cooling->table_solar_abundances[6]; /* Mg */
  element_abundance_solar[5] =
      Z[7] / cooling->table_solar_abundances[7]; /* Si */
  element_abundance_solar[6] =
      Z[7] / cooling->table_solar_abundances[7]; /* S  */
  element_abundance_solar[7] =
      Z[7] / cooling->table_solar_abundances[7]; /* Ca */
  element_abundance_solar[8] =
      Z[8] / cooling->table_solar_abundances[10]; /* Fe */

  // MATTHIEU: to do: Create array with pre-computed inverse abundances
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
                           const double u_cgs, const double n_H_cgs,
                           const float He_frac) {

#ifdef SWIFT_DEBUG_CHECKS
  if (u_cgs == 0.) error("Internal energy is 0. Aborting");
  if (n_H_cgs == 0.) error("Hydrogen number density is 0. Aborting");
#endif

  /* Start by recovering indices along the axis of the table */
  int index_He, index_nH, index_u;
  float delta_He_table, delta_nH_table, delta_u_table;

  /* Helium fraction axis */
  find_1d_index(cooling->table_He_frac, eagle_cooling_N_He_frac, He_frac,
                &index_He, &delta_He_table);

  /* Hydrogen number density axis */
  find_1d_index(cooling->table_nH, eagle_cooling_N_density,
                (float)log10(n_H_cgs), &index_nH, &delta_nH_table);

  /* Internal energy axis */
  find_1d_index(cooling->table_u, eagle_cooling_N_temperature,
                (float)log10(u_cgs), &index_u, &delta_u_table);

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
    return T * u_cgs / pow(10., (double)cooling->table_u[0]);
  else
    return T;
}

/**
 * @brief Compute the Compton cooling rate from the CMB at a given
 * redshift, electron abundance, temperature and Hydrogen density.
 *
 * Uses an analytic formula.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param redshift The current redshift.
 * @param n_H_cgs The Hydrogen number density in CGS units.
 * @param temperature The temperature.
 * @param electron_abundance The electron abundance.
 */
double eagle_Compton_cooling_rate(const struct cooling_function_data* cooling,
                                  const double redshift, const double n_H_cgs,
                                  const double temperature,
                                  const double electron_abundance) {

  const double zp1 = 1. + redshift;
  const double zp1p2 = zp1 * zp1;
  const double zp1p4 = zp1p2 * zp1p2;

  /* CMB temperature at this redshift */
  const double T_CMB = cooling->T_CMB_0 * zp1;

  /* Compton cooling rate */
  return cooling->compton_rate_cgs * (temperature - T_CMB) * zp1p4 *
         electron_abundance / n_H_cgs;
}

/**
 * @brief Computes the cooling rate corresponding to a given internal energy,
 * hydrogen number density, Helium fraction, redshift and metallicity from
 * all the possible channels.
 *
 * 1) Metal-free cooling:
 * We interpolate the flattened 4D table 'H_and_He_net_heating' that is
 * arranged in the following way:
 * - 1st dim: redshift, length = 2
 * - 2nd dim: Helium fraction, length = eagle_cooling_N_He_frac
 * - 3rd dim: Hydrogen density, length = eagle_cooling_N_density
 * - 4th dim: Internal energy, length = eagle_cooling_N_temperature
 *
 * 2) Electron abundance
 * We compute the electron abundance by interpolating the flattened 4d table
 * 'H_and_He_electron_abundance' that is arranged in the following way:
 * - 1st dim: redshift, length = 2
 * - 2nd dim: Helium fraction, length = eagle_cooling_N_He_frac
 * - 3rd dim: Hydrogen density, length = eagle_cooling_N_density
 * - 4th dim: Internal energy, length = eagle_cooling_N_temperature
 *
 * 3) Compton cooling is applied via the analytic formula.
 *
 * 4) Solar electron abudance
 * We compute the solar electron abundance by interpolating the flattened 3d
 * table
 * 'solar_electron_abundance' that is arranged in the following way:
 * - 1st dim: redshift, length = 2
 * - 2nd dim: Hydrogen density, length = eagle_cooling_N_density
 * - 3rd dim: Internal energy, length = eagle_cooling_N_temperature
 *
 * 5) Metal-line cooling
 * For each tracked element we interpolate the flattened 4D table
 * 'table_metals_net_heating' that is arrange in the following way:
 * - 1st dim: redshift, length = 2
 * - 2nd dim: element, length = eagle_cooling_N_metal
 * - 3rd dim: Hydrogen density, length = eagle_cooling_N_density
 * - 4th dim: Internal energy, length = eagle_cooling_N_temperature
 *
 * Note that this is a fake 4D interpolation as we do not interpolate
 * along the 2nd dimension. We just do this once per element.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param u_cgs The internal energy in CGS units.
 * @param n_H_cgs The Hydrogen number density in CGS units.
 * @param He_frac The Helium fraction.
 * @param redshift The current redshift.
 * @param Z The element metallicity expressed as element mass fractions.
 */
double eagle_total_cooling_rate(const struct cooling_function_data* cooling,
                                const double u_cgs, const double n_H_cgs,
                                const float He_frac, const double redshift,
                                const float Z[chemistry_element_count]) {

#ifdef SWIFT_DEBUG_CHECKS
  if (u_cgs == 0.) error("Internal energy is 0. Aborting");
  if (n_H_cgs == 0.) error("Hydrogen number density is 0. Aborting");
#endif

  /* Get the temperature corresponding to this energy */
  const float temperature =
      eagle_convert_u_to_T(cooling, u_cgs, n_H_cgs, He_frac);

  /* Start by recovering indices along the axis of the tables */
  int index_T, index_He, index_nH;
  float delta_T_table, delta_He_table, delta_nH_table;

  /* Temperature axis */
  find_1d_index(cooling->table_temperatures, eagle_cooling_N_temperature,
                log10f(temperature), &index_T, &delta_T_table);

  /* Helium fraction axis */
  find_1d_index(cooling->table_He_frac, eagle_cooling_N_He_frac, He_frac,
                &index_He, &delta_He_table);

  /* Hydrogen number density axis */
  find_1d_index(cooling->table_nH, eagle_cooling_N_density,
                (float)log10(n_H_cgs), &index_nH, &delta_nH_table);

  /* For the redshift, we can use the pre-computed delta stored in the cooling
   * function structure */
  float delta_z_table = cooling->delta_z_table;

  /**********************/
  /* Electron abundance */
  /**********************/

  const double electron_abundance = interpolation_4d(
      cooling->table_H_and_He_electron_abundance, 0, index_He, index_nH,
      index_T, 2, eagle_cooling_N_He_frac, eagle_cooling_N_density,
      eagle_cooling_N_temperature, delta_z_table, delta_He_table,
      delta_nH_table, delta_T_table);

  /**********************/
  /* Metal-free cooling */
  /**********************/

  double Lambda_net = interpolation_4d(
      cooling->table_H_and_He_net_heating, 0, index_He, index_nH, index_T, 2,
      eagle_cooling_N_He_frac, eagle_cooling_N_density,
      eagle_cooling_N_temperature, delta_z_table, delta_He_table,
      delta_nH_table, delta_T_table);

#ifdef SWIFT_DEBUG_CHECKS
  if (Lambda_net == 0.) error("Cooling is 0");
#endif

  /**********************/
  /* Compton cooling    */
  /**********************/

  /* Do we need to add the inverse Compton cooling? */
  /* It is *not* stored in the tables before re-ionisation */
  if ((redshift > cooling->table_redshifts[eagle_cooling_N_redshifts - 1]) ||
      (redshift > cooling->H_reion_z)) {

    Lambda_net -= eagle_Compton_cooling_rate(cooling, redshift, n_H_cgs,
                                             temperature, electron_abundance);
  }

  /*******************************/
  /* Solar electron abundance    */
  /*******************************/
  /* Compute the electron abundances for the Solar values */
  const double solar_electron_abundance = interpolation_3d(
      cooling->table_solar_electron_abundance, 0, index_nH, index_T, 2,
      eagle_cooling_N_density, eagle_cooling_N_temperature, delta_z_table,
      delta_nH_table, delta_T_table);

  const double abundance_ratio = electron_abundance / solar_electron_abundance;

  /**********************/
  /* Metal-line cooling */
  /**********************/

  /* We start by setting the solar abundances */
  float element_abundance_solar[eagle_cooling_N_metal];
  eagle_set_cooling_abundances(cooling, Z, element_abundance_solar);

  /* Loop over the metals */
  for (int i = 0; i < eagle_cooling_N_metal; ++i) {
    if (element_abundance_solar[i] > 0.) {

      // MATTHIEU: to do: Transpose array to have Z as 1st dimension

      const double Lambda_elem = interpolation_4d_no_y(
          cooling->table_metals_net_heating, 0, i, index_nH, index_T, 2,
          eagle_cooling_N_metal, eagle_cooling_N_density,
          eagle_cooling_N_temperature, delta_z_table, 0.f, delta_nH_table,
          delta_T_table);

      Lambda_net += Lambda_elem * abundance_ratio * element_abundance_solar[i];
    }
  }

  /* Ok, we are done ! */
  return Lambda_net;
}

/**
 * @brief Compute the new energy of a particle based on its current
 * properties.
 *
 * We want to compute u_new such that u_new = u_old + dt * du/dt(u_new, X),
 * where X stands for the metallicity, density and redshift. These are
 * kept constant.
 *
 * We first compute du/dt(u_old). If dt * du/dt(u_old) is small enough, we
 * use an explicit integration and use this as our solution.
 *
 * Otherwise, we try to find a solution to the implicit time-integration
 * problem. This leads to the root-finding problem:
 *
 * f(u_new) = u_new - u_old - dt * du/dt(u_new, X) = 0
 *
 * This is done by first bracketing the solution and then iterating
 * towards the solution by reducing the window down to a certain tolerance.
 * Note there is always at least one solution since
 * f(+inf) is < 0 and f(-inf) is > 0.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param u_old_cgs The old internal energy per unit mass in CGS units.
 * @param rho_cgs The current density in CGS units.
 * @param dt_cgs The time-step over which to cool in CGS units.
 * @param delta_z The time-step expressed as a change in redshift. This is
 * negative.
 * @param redshift The current redshift.
 * @param Z The element metallicity expressed as element mass fractions.
 *
 * @return The new internal energy in CGS units.
 */
double eagle_do_cooling(const struct cooling_function_data* cooling,
                        const double u_old_cgs, const double rho_cgs,
                        const double dt_cgs, const double delta_z,
                        const double redshift,
                        const float Z[chemistry_element_count]) {

#ifdef SWIFT_DEBUG_CHECKS
  if (delta_z > 0.f) error("Invalid value for delta_z. Should be negative.");
  if (rho_cgs == 0.) error("Density is 0. Aborting.");
  if (u_old_cgs == 0.) error("Internal energy is 0. Aborting.");
#endif

  /* Hydrogen fraction */
  const float XH = Z[chemistry_element_H];

#ifdef SWIFT_DEBUG_CHECKS
  if (XH == 0.f) error("Hydrogen fraction is 0. Something is really odd...");
  if (XH > 1.f) error("Hydrogen fraction is >1. Something is really odd...");
#endif

  /* Helium fraction */
  const float He_frac =
      Z[chemistry_element_He] / (XH + Z[chemistry_element_He]);

#ifdef SWIFT_DEBUG_CHECKS
  if (He_frac > 1.f) error("Helium fraction is >1. Something is really odd...");
#endif

  /* Hydrogen number density */
  const double n_H_cgs = rho_cgs * XH / cooling->const_proton_mass_cgs;

  /* Cooling rate factor n_H^2 / rho = n_H * XH / m_p */
  const double rate_factor_cgs = n_H_cgs * XH / cooling->const_proton_mass_cgs;

  /* Compute the extra heat injected by Helium re-ionization */
  const double u_reion_cgs =
      eagle_cooling_helium_reion_extra_heat(redshift, delta_z, cooling);

  /* Turn this into a cooling (heating) rate */
  const double Lambda_reion = u_reion_cgs / (dt_cgs * rate_factor_cgs);

  /*************************************/
  /* Let's get a first guess           */
  /*************************************/

  /* First guess */
  double Lambda_net =
      Lambda_reion + eagle_total_cooling_rate(cooling, u_old_cgs, n_H_cgs,
                                              He_frac, redshift, Z);

  const double first_delta_u_cgs = rate_factor_cgs * Lambda_net * dt_cgs;

  /* Small cooling rate: Use explicit solution and abort */
  if (fabs(first_delta_u_cgs) < eagle_cooling_explicit_tolerance * u_old_cgs) {

    /* Explicit solution */
    return u_old_cgs + first_delta_u_cgs;
  }

  message("Implicit!");

  // MATTHIEU: to do: Bring forward the common calculation of indices and metal
  // arrays.

  /*************************************/
  /* Let's try to bracket the solution */
  /*************************************/

  double u_lower_cgs = u_old_cgs;
  double u_upper_cgs = u_old_cgs;

  if (first_delta_u_cgs > 0.) { /* Net heating case */

    u_upper_cgs *= sqrt(eagle_cooling_bracketing_factor);
    u_lower_cgs /= sqrt(eagle_cooling_bracketing_factor);

    /* Compute a new rate */
    double Lambda_new = eagle_total_cooling_rate(cooling, u_upper_cgs, n_H_cgs,
                                                 He_frac, redshift, Z);
    double delta_new_cgs = u_reion_cgs + Lambda_new * rate_factor_cgs * dt_cgs;

    while (u_upper_cgs - u_old_cgs - delta_new_cgs < 0.) {

      u_upper_cgs *= eagle_cooling_bracketing_factor;
      u_lower_cgs *= eagle_cooling_bracketing_factor;

      /* Compute a new rate */
      Lambda_new = eagle_total_cooling_rate(cooling, u_upper_cgs, n_H_cgs,
                                            He_frac, redshift, Z);
      delta_new_cgs = u_reion_cgs + Lambda_new * rate_factor_cgs * dt_cgs;
    }
  }

  else if (first_delta_u_cgs < 0.) { /* Net cooling case */

    message("net cooling");

    u_upper_cgs *= sqrt(eagle_cooling_bracketing_factor);
    u_lower_cgs /= sqrt(eagle_cooling_bracketing_factor);

    /* Compute a new rate */
    double Lambda_new = eagle_total_cooling_rate(cooling, u_lower_cgs, n_H_cgs,
                                                 He_frac, redshift, Z);
    double delta_new_cgs = u_reion_cgs + Lambda_new * rate_factor_cgs * dt_cgs;

    while (u_lower_cgs - u_old_cgs - delta_new_cgs > 0.) {

      u_upper_cgs /= eagle_cooling_bracketing_factor;
      u_lower_cgs /= eagle_cooling_bracketing_factor;

      /* Compute a new rate */
      Lambda_new = eagle_total_cooling_rate(cooling, u_lower_cgs, n_H_cgs,
                                            He_frac, redshift, Z);
      delta_new_cgs = u_reion_cgs + Lambda_new * rate_factor_cgs * dt_cgs;
    }
  }

  message("upper=%e lower=%e", u_upper_cgs, u_lower_cgs);

  /********************************************/
  /* We now have an upper and lower bound.    */
  /* Let's iterate by reducing the bracketing */
  /********************************************/

  double delta_upper_lower = 0.;
  double u_new_cgs;

  do {

    /* New guess */
    u_new_cgs = 0.5f * (u_lower_cgs + u_upper_cgs);

    /* New rate */
    Lambda_net =
        Lambda_reion + eagle_total_cooling_rate(cooling, u_new_cgs, n_H_cgs,
                                                He_frac, redshift, Z);

    /* New change in energy */
    const double delta_u_cgs = rate_factor_cgs * Lambda_net * dt_cgs;

    /* Bracketing */
    if (u_new_cgs - u_old_cgs - delta_u_cgs > 0.)
      u_upper_cgs = u_new_cgs;
    else
      u_lower_cgs = u_new_cgs;

    /* Width of the range */
    delta_upper_lower = u_upper_cgs - u_lower_cgs;

  } while (fabs(delta_upper_lower / u_new_cgs) >
           eagle_cooling_implicit_tolerance);

  return u_new_cgs;
}
