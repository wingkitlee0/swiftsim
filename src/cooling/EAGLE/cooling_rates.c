
/* Local headers */
#include "cooling.h"
#include "cooling_interpolation.h"
#include "cooling_tables.h"

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

float eagle_do_cooling(const struct cooling_function_data* cooling,
                       const float uold_cgs, const float rho_cgs,
                       const float dt_cgs, const float delta_z,
                       const float redshift,
                       const float Z[chemistry_element_count]) {

  /* Hydrogen fraction */
  const float XH = Z[chemistry_element_H];

  /* Compute Helium fraction */
  const float He_frac =
      Z[chemistry_element_He] / (XH + Z[chemistry_element_He]);

  /* Hydrogen number density */
  const float n_H_cgs = rho_cgs * XH / cooling->const_proton_mass_cgs;

  /* Cooling rate factor n_H^2 / rho = n_H * XH / m_p */
  const float rate_factor_cgs = n_H_cgs * XH / cooling->const_proton_mass_cgs;

  /* Compute the extra heat injected by Helium re-ionization */
  const float u_reion =
      eagle_cooling_helium_reion_extra_heat(redshift, delta_z, cooling);
  const float Lambda_reion = u_reion / (dt_cgs * rate_factor_cgs);

  return Lambda_reion + He_frac;  // MATTHIEU: temp place-holder
}
