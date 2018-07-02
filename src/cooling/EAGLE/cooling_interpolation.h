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
#ifndef SWIFT_COOLING_INTERPOLATION_EAGLE_H
#define SWIFT_COOLING_INTERPOLATION_EAGLE_H

/*! Number of different bins along the redhsift axis of the tables */
#define eagle_cooling_N_redshifts 49

/*! Number of different bins along the temperature axis of the tables */
#define eagle_cooling_N_temperature 176

/*! Number of different bins along the density axis of the tables */
#define eagle_cooling_N_density 41

/*! Number of different bins along the metal axis of the tables */
#define eagle_cooling_N_metal 9

/*! Number of different bins along the metal axis of the tables */
#define eagle_cooling_N_He_frac 7

/*! Number of different bins along the abundances axis of the tables */
#define eagle_cooling_N_abundances 11

/**
 * @brief Returns the 1d index of element with 2d indices x,y
 * from a flattened 2d array in row major order
 *
 * @param x, y Indices of element of interest
 * @param Nx, Ny Sizes of array dimensions
 */
__attribute__((always_inline)) INLINE int row_major_index_2d(int x, int y,
                                                             int Nx, int Ny) {
#ifdef SWIFT_DEBUG_CHECKS
  assert(x < Nx);
  assert(y < Ny);
#endif
  return x * Ny + y;
}

/**
 * @brief Returns the 1d index of element with 3d indices x,y,z
 * from a flattened 3d array in row major order
 *
 * @param x, y, z Indices of element of interest
 * @param Nx, Ny, Nz Sizes of array dimensions
 */
__attribute__((always_inline)) INLINE int row_major_index_3d(int x, int y,
                                                             int z, int Nx,
                                                             int Ny, int Nz) {
#ifdef SWIFT_DEBUG_CHECKS
  assert(x < Nx);
  assert(y < Ny);
  assert(z < Nz);
#endif
  return x * Ny * Nz + y * Nz + z;
}

/**
 * @brief Returns the 1d index of element with 4d indices x,y,z,w
 * from a flattened 4d array in row major order
 *
 * @param x, y, z, w Indices of element of interest
 * @param Nx, Ny, Nz, Nw Sizes of array dimensions
 */
__attribute__((always_inline)) INLINE int row_major_index_4d(int x, int y,
                                                             int z, int w,
                                                             int Nx, int Ny,
                                                             int Nz, int Nw) {
#ifdef SWIFT_DEBUG_CHECKS
  assert(x < Nx);
  assert(y < Ny);
  assert(z < Nz);
  assert(w < Nw);
#endif
  return x * Ny * Nz * Nw + y * Nz * Nw + z * Nw + w;
}

/**
 * @brief Finds the index of a value in a table and compute delta to nearest
 * element.
 *
 * This function assumes the table is monotonically increasing with a constant
 * difference between adjacent values.
 *
 * The returned difference is expressed in units of the table separation. This
 * means dx = (x - table[i]) / (table[i+1] - table[i]). It is always between
 * 0 and 1.
 *
 * @param table The table to search in.
 * @param size The number of elements in the table.
 * @param x The value to search for.
 * @param i (return) The index in the table of the element.
 * @param *dx (return) The difference between x and table[i]
 */
__attribute__((always_inline)) INLINE void find_1d_index(const float *table,
                                                         const int size,
                                                         const float x, int *i,
                                                         float *dx) {

  const float delta = (size - 1) / (table[size - 1] - table[0]);

  // MATTHIEU: to do: Exploit alignment of the arrays

  if (x < table[0]) { /* We are below the first element */
    *i = 0;
    *dx = 0.f;
  } else if (x >= table[size - 1]) { /* We are after the last element */
    *i = size - 2;
    *dx = 1.f;
  } else { /* Normal case */
    *i = (x - table[0]) * delta;
    *dx = (x - table[*i]) * delta;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (*dx < -0.001f || *dx > 1.001f) error("Invalid distance found dx=%e", *dx);
#endif
}

/**
 * @brief Interpolate a flattened 3D table at a given position.
 *
 * This function uses linear interpolation along each axis.
 *
 * @param table The 3D table to interpolate.
 * @param xi, yi, zi Indices of element of interest.
 * @param Nx, Ny, Nz Sizes of array dimensions.
 * @param dx, dy, dz Distance between the point and the index in units of
 * the grid spacing.
 */
__attribute__((always_inline)) INLINE float interpolation_3d(
    const float *table, const int xi, const int yi, const int zi, const int Nx,
    const int Ny, const int Nz, const float dx, const float dy,
    const float dz) {

#ifdef SWIFT_DEBUG_CHECKS
  if (dx < -0.001f || dx > 1.001f) error("Invalid dx=%e", dx);
  if (dy < -0.001f || dy > 1.001f) error("Invalid dx=%e", dy);
  if (dz < -0.001f || dz > 1.001f) error("Invalid dx=%e", dz);
#endif

  const float tx = 1.f - dx;
  const float ty = 1.f - dy;
  const float tz = 1.f - dz;

  // MATTHIEU: to do: re-arrange to exploit faster access on the last entries
  // MATTHIEU: to do: Exploit alignment of the arrays

  /* Linear interpolation along each axis. We read the table 2^3=8 times */
  float result = tx * ty * tz *
                 table[row_major_index_3d(xi + 0, yi + 0, zi + 0, Nx, Ny, Nz)];

  result += tx * ty * dz *
            table[row_major_index_3d(xi + 0, yi + 0, zi + 1, Nx, Ny, Nz)];
  result += tx * dy * tz *
            table[row_major_index_3d(xi + 0, yi + 1, zi + 0, Nx, Ny, Nz)];
  result += dx * ty * tz *
            table[row_major_index_3d(xi + 1, yi + 0, zi + 0, Nx, Ny, Nz)];

  result += tx * dy * dz *
            table[row_major_index_3d(xi + 0, yi + 1, zi + 1, Nx, Ny, Nz)];
  result += dx * ty * dz *
            table[row_major_index_3d(xi + 1, yi + 0, zi + 1, Nx, Ny, Nz)];
  result += dx * dy * tz *
            table[row_major_index_3d(xi + 1, yi + 1, zi + 0, Nx, Ny, Nz)];

  result += dx * dy * dz *
            table[row_major_index_3d(xi + 1, yi + 1, zi + 1, Nx, Ny, Nz)];

  return result;
}

/**
 * @brief Interpolate a flattened 4D table at a given position.
 *
 * This function uses linear interpolation along each axis.
 *
 * @param table The 4D table to interpolate.
 * @param xi, yi, zi, wi Indices of element of interest.
 * @param Nx, Ny, Nz, Nw Sizes of array dimensions.
 * @param dx, dy, dz, dw Distance between the point and the index in units of
 * the grid spacing.
 */
__attribute__((always_inline)) INLINE float interpolation_4d(
    const float *table, const int xi, const int yi, const int zi, const int wi,
    const int Nx, const int Ny, const int Nz, const int Nw, const float dx,
    const float dy, const float dz, const float dw) {

#ifdef SWIFT_DEBUG_CHECKS
  if (dx < -0.001f || dx > 1.001f) error("Invalid dx=%e", dx);
  if (dy < -0.001f || dy > 1.001f) error("Invalid dx=%e", dy);
  if (dz < -0.001f || dz > 1.001f) error("Invalid dx=%e", dz);
  if (dw < -0.001f || dw > 1.001f) error("Invalid dx=%e", dw);
#endif

  const float tx = 1.f - dx;
  const float ty = 1.f - dy;
  const float tz = 1.f - dz;
  const float tw = 1.f - dw;

  // MATTHIEU: to do: re-arrange to exploit faster access on the last entries
  // MATTHIEU: to do: Exploit alignment of the arrays

  /* Linear interpolation along each axis. We read the table 2^4=16 times */
  float result =
      tx * ty * tz * tw *
      table[row_major_index_4d(xi + 0, yi + 0, zi + 0, wi + 0, Nx, Ny, Nz, Nw)];

  result +=
      tx * ty * tz * dw *
      table[row_major_index_4d(xi + 0, yi + 0, zi + 0, wi + 1, Nx, Ny, Nz, Nw)];
  result +=
      tx * ty * dz * tw *
      table[row_major_index_4d(xi + 0, yi + 0, zi + 1, wi + 0, Nx, Ny, Nz, Nw)];
  result +=
      tx * dy * tz * tw *
      table[row_major_index_4d(xi + 0, yi + 1, zi + 0, wi + 0, Nx, Ny, Nz, Nw)];
  result +=
      dx * ty * tz * tw *
      table[row_major_index_4d(xi + 1, yi + 0, zi + 0, wi + 0, Nx, Ny, Nz, Nw)];

  result +=
      tx * ty * dz * dw *
      table[row_major_index_4d(xi + 0, yi + 0, zi + 1, wi + 1, Nx, Ny, Nz, Nw)];
  result +=
      tx * dy * tz * dw *
      table[row_major_index_4d(xi + 0, yi + 1, zi + 0, wi + 1, Nx, Ny, Nz, Nw)];
  result +=
      dx * ty * tz * dw *
      table[row_major_index_4d(xi + 1, yi + 0, zi + 0, wi + 1, Nx, Ny, Nz, Nw)];
  result +=
      tx * dy * dz * tw *
      table[row_major_index_4d(xi + 0, yi + 1, zi + 1, wi + 0, Nx, Ny, Nz, Nw)];
  result +=
      dx * ty * dz * tw *
      table[row_major_index_4d(xi + 1, yi + 0, zi + 1, wi + 0, Nx, Ny, Nz, Nw)];
  result +=
      dx * dy * tz * tw *
      table[row_major_index_4d(xi + 1, yi + 1, zi + 0, wi + 0, Nx, Ny, Nz, Nw)];

  result +=
      dx * dy * dz * tw *
      table[row_major_index_4d(xi + 1, yi + 1, zi + 1, wi + 0, Nx, Ny, Nz, Nw)];
  result +=
      dx * dy * tz * dw *
      table[row_major_index_4d(xi + 1, yi + 1, zi + 0, wi + 1, Nx, Ny, Nz, Nw)];
  result +=
      dx * ty * dz * dw *
      table[row_major_index_4d(xi + 1, yi + 0, zi + 1, wi + 1, Nx, Ny, Nz, Nw)];
  result +=
      tx * dy * dz * dw *
      table[row_major_index_4d(xi + 0, yi + 1, zi + 1, wi + 1, Nx, Ny, Nz, Nw)];

  result +=
      dx * dy * dz * dw *
      table[row_major_index_4d(xi + 1, yi + 1, zi + 1, wi + 1, Nx, Ny, Nz, Nw)];

  return result;
}

#endif /* SWIFT_COOLING_INTERPOLATION_EAGLE_H */
