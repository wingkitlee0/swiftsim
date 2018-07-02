
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

float eagle_cooling_rate(const struct cooling_function_data* cooling,
                         const float u_cgs, const float n_H_cgs,
                         const float redshift, const float delta_z,
                         const float Z[chemistry_element_count]) {

  return 1.f;
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
      Lambda_reion +
      eagle_cooling_rate(cooling, uold_cgs, n_H_cgs, redshift, delta_z, Z);

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
    float Lambda_new =
        eagle_cooling_rate(cooling, u_upper_cgs, n_H_cgs, redshift, delta_z, Z);
    float delta_new_cgs = u_reion_cgs + Lambda_new * rate_factor_cgs * dt_cgs;

    while (u_upper_cgs - uold_cgs - delta_new_cgs < 0.f) {

      u_upper_cgs *= eagle_cooling_bracketing_factor;
      u_lower_cgs *= eagle_cooling_bracketing_factor;

      /* Compute a new rate */
      Lambda_new = eagle_cooling_rate(cooling, u_upper_cgs, n_H_cgs, redshift,
                                      delta_z, Z);
      delta_new_cgs = u_reion_cgs + Lambda_new * rate_factor_cgs * dt_cgs;
    }
  }

  if (delta_u_cgs < 0.f) { /* Net cooling case */

    u_upper_cgs *= sqrtf(eagle_cooling_bracketing_factor);
    u_lower_cgs /= sqrtf(eagle_cooling_bracketing_factor);

    /* Compute a new rate */
    float Lambda_new =
        eagle_cooling_rate(cooling, u_lower_cgs, n_H_cgs, redshift, delta_z, Z);
    float delta_new_cgs = u_reion_cgs + Lambda_new * rate_factor_cgs * dt_cgs;

    while (u_lower_cgs - uold_cgs - delta_new_cgs > 0.f) {

      u_upper_cgs /= eagle_cooling_bracketing_factor;
      u_lower_cgs /= eagle_cooling_bracketing_factor;

      /* Compute a new rate */
      Lambda_new = eagle_cooling_rate(cooling, u_lower_cgs, n_H_cgs, redshift,
                                      delta_z, Z);
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
    Lambda_net = Lambda_reion + eagle_cooling_rate(cooling, u_cgs, n_H_cgs,
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
