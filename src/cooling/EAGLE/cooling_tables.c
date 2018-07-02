/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#include "config.h"

/* Local headers */
#include "cooling.h"
#include "cooling_interpolation.h"
#include "cooling_tables.h"

/* Library headers */
#include <hdf5.h>

/*! Name of the elements in the order they are storeed in the files */
static const char* eagle_tables_element_name[eagle_cooling_N_metal] = {
    "Carbon",  "Nitrogen", "Oxygen",  "Neon", "Magnesium",
    "Silicon", "Sulphur",  "Calcium", "Iron"};

/**
 * @brief Returns the index of the current redshift in the cooling table.
 * Also compute the delta redshift since the last table.
 *
 * @param z The current redshift.
 * @param cooling The #cooling_function_data used in the run.
 * @param z_index (return) The index along the redshift axis of the cooling
 * table.
 * @param delta_z (return) The difference in redshift to the table z_index.
 */
void eagle_get_redshift_table_index(const float z,
                                    const struct cooling_function_data* cooling,
                                    int* z_index, float* delta_z) {

  const float H_reion_z = cooling->H_reion_z;
  const float* table_redshifts = cooling->table_redshifts;

  /* Before the redshift of reionisation */
  if (z > H_reion_z) {

    *z_index = eagle_cooling_N_redshifts;
    *delta_z = 0.f;
  }

  /* Between reionisation and first table */
  else if ((z <= H_reion_z) &&
           (z > table_redshifts[eagle_cooling_N_redshifts - 1])) {

    *z_index = eagle_cooling_N_redshifts + 1;
    *delta_z = 0.f;
  }

  /* After the last table */
  else if (z < table_redshifts[0]) {

    *z_index = 0;
    *delta_z = 0.f;
  }

  /* Normal case: figure out where we are */
  else {

    /* Search the list of redshifts */
    for (int i = eagle_cooling_N_redshifts - 2; i >= 0; --i) {

      if (z > table_redshifts[i]) {

        *z_index = i;
        *delta_z = (z - table_redshifts[i]) /
                   (table_redshifts[i + 1] - table_redshifts[i]);

        break;
      }
    }
  }
}

/**
 * @brief Read the redshifts of the cooling tables from a file.
 *
 * Also checks that the redshift are in increasing order.
 *
 * @param cooling The cooling data we play with.
 */
void eagle_cooling_init_redshift_tables(struct cooling_function_data* cooling) {

  /* Allocate memory */
  cooling->table_redshifts = malloc(eagle_cooling_N_redshifts * sizeof(float));
  if (cooling->table_redshifts == NULL)
    error("Error allocating table of redshifts");

  // MATTHIEU: Change this to read a generic file name

  /* Open file */
  FILE* file = fopen("redshifts.dat", "r");
  if (file == NULL) error("Impossible to open the list of redshifts");

  int num_redshifts;
  if (fscanf(file, "%d", &num_redshifts) != 1)
    error("Error reading the number of redshifts from the file");

  /* Check that the number of redshifts matches the code. */
  if (num_redshifts != eagle_cooling_N_redshifts)
    error("Number of redshifts in the file does not match the code.");

  /* Read the redshifts */
  for (int i = 0; i < eagle_cooling_N_redshifts; ++i) {

    float buffer;
    if (fscanf(file, "%f", &buffer) != 1)
      error("Impossible to read redshift number %d", i);

    cooling->table_redshifts[i] = buffer;
  }

  /* Check that the redshifts are in increasing order */
  for (int i = 0; i < eagle_cooling_N_redshifts - 1; ++i)
    if (cooling->table_redshifts[i] > cooling->table_redshifts[i + 1])
      error("Redhsift tables are not in increasing order");

  fclose(file);
}

/**
 * @brief Read an HDF5 integer
 *
 * This function make no check that the types match.
 *
 * @param h_file The (opened) HDF5 file.
 * @param name The name of the integer to read.
 */
int read_hdf5_int(const hid_t h_file, const char* name) {

  int value;

  const hid_t h_set = H5Dopen(h_file, name, H5P_DEFAULT);
  if (h_set < 0) error("Error while opening data space for '%s'", name);
  const hid_t h_err =
      H5Dread(h_set, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
  if (h_err < 0) error("Error while reading data space for '%s'", name);
  H5Dclose(h_set);

  return value;
}

/**
 * @brief Read an HDF5 array into an allocated array.
 *
 * This function make no check that the array sizes or types match.
 *
 * @param h_file The (opened) HDF5 file.
 * @param name The name of the array to read.
 * @param array The array to fill.
 */
void read_hdf5_array(const hid_t h_file, const char* name, float* array) {

  const hid_t h_set = H5Dopen(h_file, name, H5P_DEFAULT);
  if (h_set < 0) error("Error while opening data space for '%s'", name);
  const hid_t h_err =
      H5Dread(h_set, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, array);
  if (h_err < 0) error("Error while reading data space for '%s'", name);
  H5Dclose(h_set);
}

/**
 * @brief Read the header of the cooling table file.
 *
 * This reads all the arrays that are independant of redshift and allocates
 * all the arrays.
 *
 * @param cooling The #cooling_function_data we play with.
 */
void eagle_read_cooling_table_header(struct cooling_function_data* cooling) {

  // MATTHIEU: Change this to read a generic file name
  const char* filename = "z_0.000.hdf5";

  /* Open file */
  const hid_t h_file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (h_file < 0) error("Error while opening file '%s'.", filename);

  /* Start by checking that the table dimensions match */
  int temp;

  /* Number of density bins */
  temp = read_hdf5_int(h_file, "/Header/Number_of_density_bins");
  if (temp != eagle_cooling_N_density)
    error("Number of density bins in the file does not match the code.");
  message("Cooling table has %d elements along the density axis", temp);

  /* Number of temperature bins */
  temp = read_hdf5_int(h_file, "/Header/Number_of_temperature_bins");
  if (temp != eagle_cooling_N_temperature)
    error("Number of temperature bins in the file does not match the code.");
  message("Cooling table has %d elements along the temperature axis", temp);

  /* Number of helium_fractions */
  temp = read_hdf5_int(h_file, "/Header/Number_of_helium_fractions");
  if (temp != eagle_cooling_N_He_frac)
    error("Number of helium fractions in the file does not match the code.");
  message("Cooling table has %d helium fractions", temp);

  /* Number of abundances */
  temp = read_hdf5_int(h_file, "/Header/Abundances/Number_of_abundances");
  if (temp != eagle_cooling_N_abundances)
    error("Number of abundances in the file does not match the code.");
  message("Cooling table has %d abundances", temp);

  /* Number of metals */
  temp = read_hdf5_int(h_file, "/Header/Number_of_metals");
  if (temp != eagle_cooling_N_metal)
    error("Number of metals in the file does not match the code.");
  message("Cooling table has %d metals", temp);

  /* Now allocate all the header arrays */

  // MATTHIEU: to do: Use aligned allocation.

  /* Table for log_10(T) */
  cooling->table_temperatures =
      malloc(eagle_cooling_N_temperature * sizeof(float));
  if (cooling->table_temperatures == NULL)
    error("Error allocating memory for temperatures table");

  /* Table for log_10(n_H) */
  cooling->table_nH = malloc(eagle_cooling_N_density * sizeof(float));
  if (cooling->table_nH == NULL)
    error("Error allocating memory for H density table");

  /* Table for log_10(u) */
  cooling->table_u = malloc(eagle_cooling_N_temperature * sizeof(float));
  if (cooling->table_u == NULL)
    error("Error allocating memory for internal energy table");

  /* Table for He abundance */
  cooling->table_He_frac = malloc(eagle_cooling_N_He_frac * sizeof(float));
  if (cooling->table_He_frac == NULL)
    error("Error allocating memory for Helium abundance table");

  /* Table for solar abundance */
  cooling->table_solar_abundances =
      malloc(eagle_cooling_N_abundances * sizeof(float));
  if (cooling->table_solar_abundances == NULL)
    error("Error allocating memory for solar abundance table");

  /* Now, read the content of these tables */

  /* Temperature bins */
  read_hdf5_array(h_file, "/Solar/Temperature_bins",
                  cooling->table_temperatures);

  /* Density bins */
  read_hdf5_array(h_file, "/Solar/Hydrogen_density_bins", cooling->table_nH);

  /* Internal energy bins */
  read_hdf5_array(h_file, "/Metal_free/Temperature/Energy_density_bins",
                  cooling->table_u);

  /* Helium fraction bins */
  read_hdf5_array(h_file, "/Metal_free/Helium_mass_fraction_bins",
                  cooling->table_He_frac);

  /* Solar abundances bins */
  read_hdf5_array(h_file, "/Header/Abundances/Solar_mass_fractions",
                  cooling->table_solar_abundances);

  // MATTHIEU: to do: Check that the metals are in the same order

  H5Fclose(h_file);

  /* Finally turn the quantities in their log values */

  /* Temperatures */
  for (int i = 0; i < eagle_cooling_N_temperature; ++i)
    cooling->table_temperatures[i] = log10f(cooling->table_temperatures[i]);

  /* Densities */
  for (int i = 0; i < eagle_cooling_N_density; ++i)
    cooling->table_nH[i] = log10f(cooling->table_nH[i]);

  /* Internal energies */
  for (int i = 0; i < eagle_cooling_N_temperature; ++i)
    cooling->table_u[i] = log10f(cooling->table_u[i]);
}

/**
 * @brief Compute the solar metallicity based on the HDF5 table
 *
 * This function assumes the elements are sorted in the right order
 * such that H and He have index 0 and 1 respectively.
 *
 * @param cooling The #cooling_function_data we play with.
 */
void eagle_set_solar_metallicity(struct cooling_function_data* cooling) {

  cooling->solar_metallicity = 1.f;

  /* Deduct Hydrogen abundance */
  cooling->solar_metallicity -= cooling->table_solar_abundances[0];

  /* Deduct Helium abundance */
  cooling->solar_metallicity -= cooling->table_solar_abundances[1];

  message("Solar metallicity of the tables: Z=%f", cooling->solar_metallicity);
}

/**
 * @brief Allocates memory for the cooling tables themselves.
 *
 * The 3D and 4D tables are allocated as a flattened 1D array.
 *
 * @param cooling The #cooling_function_data we play with.
 */
void eagle_allocate_cooling_tables(struct cooling_function_data* cooling) {

  /* Energy to temperature table */
  const size_t size_u_to_temp = eagle_cooling_N_He_frac *
                                eagle_cooling_N_density *
                                eagle_cooling_N_temperature * 2 * sizeof(float);

  if (posix_memalign((void**)&cooling->table_u_to_temp, 64, size_u_to_temp) !=
      0)
    error("Error allocating memory for the u->T  table");
  bzero(cooling->table_u_to_temp, size_u_to_temp);

  /* Solar electron abundance table */
  const size_t size_e_abundance =
      eagle_cooling_N_density * eagle_cooling_N_temperature * 2 * sizeof(float);

  if (posix_memalign((void**)&cooling->table_solar_electron_abundance, 64,
                     size_e_abundance) != 0)
    error("Error allocating memory for the solar electron abundance table");
  bzero(cooling->table_solar_electron_abundance, size_e_abundance);

  /* Metal cooling table */
  const size_t size_metal_cool =
      eagle_cooling_N_metal * eagle_cooling_N_density *
      eagle_cooling_N_temperature * 2 * sizeof(float);

  if (posix_memalign((void**)&cooling->table_metals_net_heating, 64,
                     size_metal_cool) != 0)
    error("Error allocating memory for the metal cooling table");
  bzero(cooling->table_metals_net_heating, size_metal_cool);

  /* H and He cooling table */
  const size_t size_HandHe_cool =
      eagle_cooling_N_He_frac * eagle_cooling_N_density *
      eagle_cooling_N_temperature * 2 * sizeof(float);

  if (posix_memalign((void**)&cooling->table_H_and_He_net_heating, 64,
                     size_HandHe_cool) != 0)
    error("Error allocating memory for the H and He cooling table");
  bzero(cooling->table_H_and_He_net_heating, size_HandHe_cool);

  /* H and He abundances */
  const size_t size_HandHe_abundances =
      eagle_cooling_N_He_frac * eagle_cooling_N_density *
      eagle_cooling_N_temperature * 2 * sizeof(float);

  if (posix_memalign((void**)&cooling->table_H_and_He_electron_abundance, 64,
                     size_HandHe_abundances) != 0)
    error("Error allocating memory for the H and He abundances table");
  bzero(cooling->table_H_and_He_electron_abundance, size_HandHe_abundances);
}

/**
 * @brief Read one set of tables from an HDF5 file.
 *
 * @param cooling The #cooling_function_data we play with.
 * @param filename The HDF5 file to read from.
 * @param z_index The index of this table along the redshift axis of the tables.
 */
void eagle_read_one_table(struct cooling_function_data* cooling,
                          const char* filename, int z_index) {

  message("Reading table '%s'", filename);

  /* Open file */
  hid_t h_file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (h_file < 0) error("Error while opening file '%s'.", filename);

  /* Temporary arrays for transactions */

  /* Metal cooling rates */
  int n_elem = eagle_cooling_N_temperature * eagle_cooling_N_density;
  float* net_cooling_rate = malloc(n_elem * sizeof(float));
  if (net_cooling_rate == NULL) error("Can't allocate temporary memory");

  /* Electron abundance */
  float* electron_abundance = malloc(n_elem * sizeof(float));
  if (electron_abundance == NULL) error("Can't allocate temporary memory");

  /* H + He cooling rate */
  n_elem = eagle_cooling_N_He_frac * eagle_cooling_N_temperature *
           eagle_cooling_N_density;
  float* he_net_cooling_rate = malloc(n_elem * sizeof(float));
  if (he_net_cooling_rate == NULL) error("Can't allocate temporary memory");

  /* u->T conversion table */
  float* temperature = malloc(n_elem * sizeof(float));
  if (temperature == NULL) error("Can't allocate temporary memory");

  /* Electron abundance table */
  float* he_electron_abundance = malloc(n_elem * sizeof(float));
  if (he_electron_abundance == NULL) error("Can't allocate temporary memory");

  /* Read the metal rates */
  for (int elem = 0; elem < eagle_cooling_N_metal; ++elem) {

    char name[200];
    sprintf(name, "/%s/Net_Cooling", eagle_tables_element_name[elem]);
    read_hdf5_array(h_file, name, net_cooling_rate);

    /* Transpose the array to make temperature the last index */
    for (int i = 0; i < eagle_cooling_N_density; ++i) {
      for (int j = 0; j < eagle_cooling_N_temperature; ++j) {

        /* Get indices in both tables */
        const int index_4d = row_major_index_4d(
            z_index, elem, i, j, 2, eagle_cooling_N_metal,
            eagle_cooling_N_density, eagle_cooling_N_temperature);

        const int index_2d = row_major_index_2d(
            j, i, eagle_cooling_N_temperature, eagle_cooling_N_density);

        /* Note the minus sign here */
        cooling->table_metals_net_heating[index_4d] =
            -net_cooling_rate[index_2d];
      }
    }
  }

  /* Read the H + He rates */
  read_hdf5_array(h_file, "/Metal_free/Net_Cooling", he_net_cooling_rate);

  /* Read the u->T conversion table */
  read_hdf5_array(h_file, "/Metal_free/Temperature/Temperature", temperature);

  /* Read the He electron abundance table */
  read_hdf5_array(h_file, "/Metal_free/Electron_density_over_n_h",
                  he_electron_abundance);

  /* Transpose the arrays to make temperature the last index */
  for (int i = 0; i < eagle_cooling_N_He_frac; ++i) {
    for (int j = 0; j < eagle_cooling_N_temperature; ++j) {
      for (int k = 0; k < eagle_cooling_N_density; ++k) {

        /* Get indices in both tables */
        const int index_4d = row_major_index_4d(
            z_index, i, k, j, 2, eagle_cooling_N_He_frac,
            eagle_cooling_N_density, eagle_cooling_N_temperature);

        const int index_3d = row_major_index_3d(
            i, j, k, eagle_cooling_N_He_frac, eagle_cooling_N_temperature,
            eagle_cooling_N_density);

        /* Note the minus sign */
        cooling->table_H_and_He_net_heating[index_4d] =
            -he_net_cooling_rate[index_3d];

        cooling->table_H_and_He_electron_abundance[index_4d] =
            he_electron_abundance[index_3d];

        cooling->table_u_to_temp[index_4d] = log10f(temperature[index_3d]);
      }
    }
  }

  /* Read the electron abundance table */
  read_hdf5_array(h_file, "/Solar/Electron_density_over_n_h",
                  electron_abundance);

  /* Transpose the arrays to make temperature the last index */
  for (int i = 0; i < eagle_cooling_N_temperature; ++i) {
    for (int j = 0; j < eagle_cooling_N_density; ++j) {

      const int index_3d =
          row_major_index_3d(z_index, j, i, 2, eagle_cooling_N_density,
                             eagle_cooling_N_temperature);

      const int index_2d = row_major_index_2d(i, j, eagle_cooling_N_temperature,
                                              eagle_cooling_N_density);

      cooling->table_solar_electron_abundance[index_3d] =
          electron_abundance[index_2d];
    }
  }

  /* Free everything */
  free(net_cooling_rate);
  free(electron_abundance);
  free(he_net_cooling_rate);
  free(temperature);
  free(he_electron_abundance);
}

/**
 * @brief Read the two cooling tables bracketing the current redshift.
 *
 * @param cooling The #cooling_function_data we play with.
 */
void eagle_load_cooling_tables(struct cooling_function_data* cooling) {

  const int low_z_index = cooling->low_z_index;
  const int high_z_index = cooling->high_z_index;

  // MATTHIEU: Optimize this. Can copy one table in the place of the other.

  char filename[200];

  /* Read the first table */
  sprintf(filename, "z_%1.3f.hdf5", cooling->table_redshifts[low_z_index]);
  eagle_read_one_table(cooling, filename, 0);

  /* Read the second table */
  sprintf(filename, "z_%1.3f.hdf5", cooling->table_redshifts[high_z_index]);
  eagle_read_one_table(cooling, filename, 1);
}

/**
 * @brief Checks the tables that are currently loaded in memory and read
 * new ones if necessary.
 *
 * @param cooling The #cooling_function_data we play with.
 * @param index_z The index along the redshift axis of the tables of the current
 * z.
 */
void eagle_check_cooling_tables(struct cooling_function_data* cooling,
                                int index_z) {

  /* Do we already have the right table in memory? */
  if (cooling->low_z_index == index_z) return;

  if (index_z >= eagle_cooling_N_redshifts) {

    error("Missing implementation");
    // MATTHIEU: Add reading of high-z tables

    /* Record the table indices */
    cooling->low_z_index = index_z;
    cooling->high_z_index = index_z;
  } else {

    /* Record the table indices */
    cooling->low_z_index = index_z;
    cooling->high_z_index = index_z + 1;

    /* Load the damn thing */
    eagle_load_cooling_tables(cooling);
  }
}
