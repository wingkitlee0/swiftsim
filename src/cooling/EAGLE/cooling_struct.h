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
#ifndef SWIFT_COOLING_STRUCT_EAGLE_H
#define SWIFT_COOLING_STRUCT_EAGLE_H

/**
 * @brief Properties of the cooling function.
 */
struct cooling_function_data {

  /* Some constants -------------------------- */

  /*! Proton mass in cgs */
  double const_proton_mass_cgs;

  /*! Temperature of the CMB at redshift 0 */
  double T_CMB_0;

  /*! Compton rate in cgs [g cm^2 s^-3 K^-1] */
  double compton_rate_cgs;

  /*! Redshift of H reionization */
  double H_reion_z;

  /*! Centre of the Gaussian used for He II reionization */
  double He_reion_z_centre;

  /*! Spread of the Gaussian used for He II reionization */
  double He_reion_z_sigma;

  /*! Energy injected per gram by He II reionization */
  double He_reion_heat_cgs;

  /* Cooling table meta-data -------------------------- */

  /*! Solar mettalicity of the tables */
  float eagle_cooling_N_abunsolar_metallicity;

  /*! The redshifts of the cooling tables */
  float *table_redshifts;

  /*! The log10 of the temperature of the cooling tables */
  float *table_temperatures;

  /*! The log10 of the H number density of the cooling tables */
  float *table_nH;

  /*! The log10 of the internal energy of the cooling tables */
  float *table_u;

  /*! The He abundances of the cooling tables */
  float *table_He_frac;

  /*! The solar abundances of the cooling tables */
  float *table_solar_abundances;

  /*! The solar metallicity of the cooling tables */
  float solar_metallicity;

  /* Cooling tables -------------------------- */

  /*! Thermal energy to temperature */
  float *table_u_to_temp;

  /*! Net heating from metals */
  float *table_metals_net_heating;

  /*! Net heating from H and He */
  float *table_H_and_He_net_heating;

  /*! Net heating from H and He */
  float *table_H_and_He_electron_abundance;

  /*! Net heating from H and He */
  float *table_solar_electron_abundance;

  /* Shared constants for this time-step --------------------------- */

  /*! The index along the redshift index of the cooling tables. */
  int index_z;

  /*! Distance in redshift from the index_z entries. */
  float delta_z_table;

  /*! Index of the low-redshift table currently loaded. */
  int low_z_index;

  /*! Index of the high-redshift table currently loaded. */
  int high_z_index;
};

/**
 * @brief Properties of the cooling stored in the extended particle data.
 */
struct cooling_xpart_data {};

#endif /* SWIFT_COOLING_STRUCT_EAGLE_H */
