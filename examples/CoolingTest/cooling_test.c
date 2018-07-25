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
  if (argc != 3) {
    message("ERROR: Need to provide redshift and abundance ratio");
    return 1;
  }

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

  /***********************************************************/
  /* Let's start by reproducing the Wiersma et al. figure 2. */
  /***********************************************************/

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

  //float Z[chemistry_element_count] = {0.752, 0.248, 0.};
  
  /* Abundance as a function of solar */
  float sum_Z = 0.f;
  for (int i = 0; i < chemistry_element_count; ++i) {
    if (i > 1) {
      Z[i] *= abundance_ratio;
    }
    sum_Z += Z[i];
  }

  /* Print abundance pattern */
  printf("[00000.1] main: Z= [ ");
  for (int i = 0; i < chemistry_element_count; ++i) {
    printf("%f ", Z[i]);
  }
  printf("] \n");
  
  /* Helium farction */
  const double He_frac = Z[1] / (Z[0] + Z[1]);
  // const double He_frac = Z[1];
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

  message("");

  /****************************************/
  /* Let's now cool a particle over time  */
  /****************************************/

  const double n_H_cgs = 1e-4;
  const double rho_cgs = n_H_cgs * func.const_proton_mass_cgs / Z[0];
  const double m_cgs = 1e6 * 1.989e33; /* ~ EAGLE particle mass */

  //const double log_u_old_cgs = 13.32;
  const double u_old_cgs = 2e13;//pow(10., log_u_old_cgs);
  const double dt_cgs = 0.75e16; /* ~300'000'000 years */
  const double delta_z =
      -0. * dt_cgs / units_cgs_conversion_factor(&us, UNIT_CONV_TIME) / cosmo.a;

  /* Get the current temperature */
  const double T_old = eagle_convert_u_to_T(&func, u_old_cgs, n_H_cgs, He_frac);

  message("n_H_cgs= %e", n_H_cgs);
  message("rho_cgs= %e", rho_cgs);
  message("dt_cgs = %e", dt_cgs);

  const double c_cgs = gas_soundspeed_from_internal_energy(rho_cgs, u_old_cgs);
  message("c_cgs= %e", c_cgs);

  const double h_cgs = 1.2348 * cbrt(m_cgs / rho_cgs);
  message("Smoothing length: h_cgs= %e", h_cgs);

  const double dt_courant_cgs = 2.f * kernel_gamma * 0.1 * h_cgs / c_cgs;
  message("dt_courant_cgs= %e", dt_courant_cgs);

  /* Cool! */
  const double u_new_cgs =
      eagle_do_cooling(&func, u_old_cgs, rho_cgs, dt_cgs, delta_z, cosmo.z, Z);

  /* Get the new temperature */
  double T_new = eagle_convert_u_to_T(&func, u_new_cgs, n_H_cgs, He_frac);

  /* Now do the same thing with an explicit solver and small dt */
  const int num_steps = 1000000;

  double u_explicit_cgs = u_old_cgs;
  const double dt_exp_cgs = dt_cgs / num_steps;
  for (int i = 0; i < num_steps; ++i) {
    u_explicit_cgs = eagle_do_cooling(&func, u_explicit_cgs, rho_cgs,
                                      dt_exp_cgs, delta_z, cosmo.z, Z);
  }

  /* Temperature from the explicit solution */
  double T_exp = eagle_convert_u_to_T(&func, u_explicit_cgs, n_H_cgs, He_frac);

  message("u_old_cgs= %e, u_new_cgs= %e u_explicit= %e", u_old_cgs, u_new_cgs,
          u_explicit_cgs);
  message("T_old= %e T_new= %e T_exp= %e", T_old, T_new, T_exp);

  const double rate = n_H_cgs * Z[0] / func.const_proton_mass_cgs;

  file = fopen("f.dat", "w");
  fprintf(file, "# u_new, f(u_new)\n");

  double u_test_cgs = u_old_cgs;
  for (int i = 0; i < num_steps; ++i) {

    u_test_cgs = u_old_cgs - 0.999 * i * u_old_cgs / num_steps;

    const double du_dt = eagle_total_cooling_rate(&func, u_test_cgs, n_H_cgs,
                                                  He_frac, cosmo.z, Z);
    const double f_u = u_test_cgs - u_old_cgs - dt_cgs * du_dt * rate;

    fprintf(file, "%e %e %e\n", u_test_cgs, f_u, du_dt * rate * dt_cgs);
  }
  fclose(file);

#endif

  message("Done.");
  return 0;
}
