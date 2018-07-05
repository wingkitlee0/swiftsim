/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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
#include "../config.h"

#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Local includes */
#include "swift.h"

int main(int argc, char* argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

/* Choke on FPEs */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  /* Check number of arguments */
  if (argc != 3) error("Need to provide redshift and abundance ratio");

  /* Read the redshift from the command line */
  double redshift;
  sscanf(argv[1], "%lf", &redshift);
  if (redshift < 0.) error("Invalid redshift value: %f", redshift);
  if (redshift == 0.) redshift = 1e-6;

  /* Read the abundance ratio from the command line */
  double abundance_ratio = 1.;
  sscanf(argv[2], "%lf", &abundance_ratio);
  if (abundance_ratio < 0.)
    error("Invalid abundance ratio value: %f", abundance_ratio);

  /* Read in the parameter file */
  struct swift_params params;
  parser_read_file("cooling_test.yml", &params);

  /* Change the redshift to the one set by the user */
  const double a = 1. / (redshift + 1.);
  char buffer[100];
  sprintf(buffer, "Cosmology:a_begin:%f", a);
  parser_set_param(&params, buffer);

  /* Initialise the system of units */
  struct unit_system us;
  units_init_from_params(&us, &params, "InternalUnitSystem");

  message("Internal unit system: U_M = %e g.", us.UnitMass_in_cgs);
  message("Internal unit system: U_L = %e cm.", us.UnitLength_in_cgs);
  message("Internal unit system: U_t = %e s.", us.UnitTime_in_cgs);
  message("Internal unit system: U_I = %e A.", us.UnitCurrent_in_cgs);
  message("Internal unit system: U_T = %e K.", us.UnitTemperature_in_cgs);
  message("");

  /* Initialise the physical constants */
  struct phys_const phys_const;
  phys_const_init(&us, &params, &phys_const);
  phys_const_print(&phys_const);
  message("");

  /* Initialise the cosmology */
  struct cosmology cosmo;
  cosmology_init(&params, &us, &phys_const, &cosmo);
  cosmology_print(&cosmo);
  message("");

  /* Initialise the cooling function */
  struct cooling_function_data func;
  cooling_init(&params, &us, &phys_const, &func);
  cooling_print(&func);
  message("");

  /* Let's updated everything to the current redshift */
  cosmology_update(&cosmo, &phys_const, 0);
  cooling_update(&phys_const, &us, &cosmo, &func);
  message("");

/* Test the EAGLE cooling */
#if defined(COOLING_EAGLE)

  /* Let's start by reproducing the Wiersma et al. figure 2. */

  /* Solar abundances */
  float Z[chemistry_element_count] = {0.70649785,     // H
                                      0.28055534,     // He
                                      0.0020665438,   // C
                                      8.3562563E-4,   // N
                                      0.0054926244,   // O
                                      0.0014144605,   // Ne
                                      5.907064E-4,    // Mg
                                      6.825874E-4,    // Si
                                      0.0011032152};  // Fe
  /* Abundance as a function of solar */
  float sum_Z = 0.f;
  for (int i = 0; i < chemistry_element_count; ++i) {
    if (i > 1) {
      Z[i] *= abundance_ratio;
    }
    sum_Z += Z[i];
  }

  /* Helium farction */
  const double He_frac = Z[1] / sum_Z;
  message("Helium fraction: He_frac= %f", He_frac);

  /* Prepare the file header */
  FILE* file = fopen("rates.txt", "w");
  fprintf(file, "# redshift = %f\n", cosmo.z);
  fprintf(file, "# Z/Z_sun = %f\n", abundance_ratio);
  fprintf(file, "# ");

  for (int j = 0; j < 7; ++j) {

    /* Hydrogen number density */
    const double log_n_H_cgs = -6. + j;
    const double n_H_cgs = pow(10., log_n_H_cgs);
    const double rho_cgs = n_H_cgs * func.const_proton_mass_cgs / Z[0];
    message("n_H = %e [cm^-3] rho = %e [g cm^-3]", n_H_cgs, rho_cgs);

    fprintf(file, "T   du/dt   ");
  }
  fprintf(file, "\n");

  /* Loop over the temperatures */
  for (int i = 0; i < 250; ++i) {

    const double log_u_cgs = 10. + i / 30.;
    const double u_cgs = pow(10., log_u_cgs);

    /* Loop over densities */
    for (int j = 0; j < 7; ++j) {

      /* Hydrogen number density */
      const double log_n_H_cgs = -6. + j;
      const double n_H_cgs = pow(10., log_n_H_cgs);

      /* Get the temperature */
      const double T = eagle_convert_u_to_T(&func, u_cgs, n_H_cgs, He_frac);

      /* Get the total cooling rate */
      const double total_du_dt_cgs =
          eagle_total_cooling_rate(&func, u_cgs, n_H_cgs, He_frac, cosmo.z, Z);

      /* Dump to file */
      fprintf(file, "%e %e ", T, total_du_dt_cgs);
    }
    fprintf(file, "\n");
  }
#endif

  return 0;
}
