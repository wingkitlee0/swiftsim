/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016   Matthieu Schaller (matthieu.schaller@durham.ac.uk).
 *               2018   Jacob Kegerreis (jacob.kegerreis@durham.ac.uk).
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
#ifndef SWIFT_TILLOTSON_EQUATION_OF_STATE_H
#define SWIFT_TILLOTSON_EQUATION_OF_STATE_H

/**
 * @file equation_of_state/planetary/tillotson.h
 *
 * Contains the Tillotson EOS functions for
 * equation_of_state/planetary/equation_of_state.h
 *
 */

/* Some standard headers. */
#include <math.h>

/* Local headers. */
#include "adiabatic_index.h"
#include "common_io.h"
#include "inline.h"
#include "units.h"
#include "physical_constants.h"
#include "equation_of_state.h"

/* ------------------------------------------------------------------------- */

// Tillotson parameters
struct Til_params {
    int mat_id;
    float rho_0, a, b, A, B, E_0, E_iv, E_cv, alpha, beta, eta_min, P_min;
    float c_TEMPORARY;
};

// Parameter values for each material (cgs units)
INLINE static void set_Til_iron(struct Til_params *mat, int mat_id) {
    mat->mat_id = mat_id;
    mat->rho_0 = 7.800;
    mat->a = 0.5;
    mat->b = 1.5;
    mat->A = 1.28e12;
    mat->B = 1.05e12;
    mat->E_0 = 9.5e10;
    mat->E_iv = 2.4e10;
    mat->E_cv = 8.67e10;
    mat->alpha = 5.0;
    mat->beta = 5.0;
    mat->eta_min = 0.0;
    mat->P_min = 0.0;

    mat->c_TEMPORARY = 9.4e-4;
}
INLINE static void set_Til_granite(struct Til_params *mat, int mat_id) {
    mat->mat_id = mat_id;
    mat->rho_0 = 2.680;
    mat->a = 0.5;
    mat->b = 1.3;
    mat->A = 1.8e11;
    mat->B = 1.8e11;
    mat->E_0 = 1.6e11;
    mat->E_iv = 3.5e10;
    mat->E_cv = 1.8e11;
    mat->alpha = 5.0;
    mat->beta = 5.0;
    mat->eta_min = 0.0;
    mat->P_min = 0.0;

    mat->c_TEMPORARY = 9.4e-4;
}
INLINE static void set_Til_water(struct Til_params *mat, int mat_id) {
    mat->mat_id = mat_id;
    mat->rho_0 = 0.998;
    mat->a = 0.7;
    mat->b = 0.15;
    mat->A = 2.18e10;
    mat->B = 1.325e11;
    mat->E_0 = 7.0e10;
    mat->E_iv = 4.19e9;
    mat->E_cv = 2.69e10;
    mat->alpha = 10.0;
    mat->beta = 5.0;
    mat->eta_min = 0.915;
    mat->P_min = 0.0;

    mat->c_TEMPORARY = 9.4e-4;
}

// Convert from cgs to internal units
INLINE static void convert_units_Til(
    struct Til_params *mat, const struct unit_system* us) {

    mat->rho_0 /= units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
    mat->A /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
    mat->B /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
    mat->E_0 /= units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
    mat->E_iv /= units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
    mat->E_cv /= units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
    mat->P_min /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);

    mat->c_TEMPORARY /= units_cgs_conversion_factor(us, UNIT_CONV_SPEED);
}

// gas_internal_energy_from_entropy
INLINE static float Til_internal_energy_from_entropy(
    float density, float entropy, struct Til_params *mat) {

    error("This EOS function is not yet implemented!");

    return 0;
}

// gas_pressure_from_entropy
INLINE static float Til_pressure_from_entropy(
    float density, float entropy, struct Til_params *mat) {

    error("This EOS function is not yet implemented!");

    return 0;
}

// gas_entropy_from_pressure
INLINE static float Til_entropy_from_pressure(
    float density, float pressure, struct Til_params *mat) {

    error("This EOS function is not yet implemented!");

    return 0;
}

// gas_soundspeed_from_entropy
INLINE static float Til_soundspeed_from_entropy(
    float density, float entropy, struct Til_params *mat) {

    error("This EOS function is not yet implemented!");

    return 0;
}

// gas_entropy_from_internal_energy
INLINE static float Til_entropy_from_internal_energy(
    float density, float u, struct Til_params *mat) {

    error("This EOS function is not yet implemented!");

    return 0;
}

// gas_pressure_from_internal_energy
INLINE static float Til_pressure_from_internal_energy(
    float density, float u, struct Til_params *mat) {

    const float eta = density / mat->rho_0;
    const float mu = eta - 1.f;
    const float nu = 1.f/eta - 1.f;
    float P_c, P_e, P;

    // Condensed or cold
    if (eta < mat->eta_min) {
        P_c = 0.f;
    }
    else {
        P_c = (mat->a + mat->b / (u / (mat->E_0 * eta*eta) + 1.f)) * density * u
            + mat->A * mu + mat->B * mu*mu;
    }
    // Expanded and hot
    P_e = mat->a*density*u + (
        mat->b * density * u / (u / (mat->E_0 * eta*eta) + 1.f)
        + mat->A*mu * exp(-mat->beta * nu)
        ) * exp(-mat->alpha * nu*nu);

    // Condensed or cold state
    if ((1.f < eta) || (u < mat->E_iv)) {
        P = P_c;
    }
    // Expanded and hot state
    else if ((eta < 1.f) && (mat->E_cv < u)) {
        P = P_e;
    }
    // Hybrid state
    else {
        P = ((u - mat->E_iv)*P_e + (mat->E_cv - u)*P_c) /
            (mat->E_cv - mat->E_iv);
    }

    // Minimum pressure
    if (P < mat->P_min) {
        P = mat->P_min;
    }

    return P;
}

// gas_internal_energy_from_pressure
INLINE static float Til_internal_energy_from_pressure(
    float density, float P, struct Til_params *mat) {

    error("This EOS function is not yet implemented!");

    return 0;
}

// gas_soundspeed_from_internal_energy
INLINE static float Til_soundspeed_from_internal_energy(
    float density, float u, struct Til_params *mat) {

//    const float eta = density / mat->rho_0;
//    const float mu = eta - 1.f;
//    const float nu = 1.f/eta - 1.f;
//    float P_c, P_e, P, c_c, c_e, c;
//
//    // Condensed or cold
//    if (eta < mat->eta_min) {
//        P_c = 0.f;
//    }
//    else {
//        P_c = (mat->a + mat->b / (u / (mat->E_0 * eta*eta) + 1.f)) * density * u
//            + mat->A * mu + mat->B * mu*mu;
//    }
//    c_c = mat->a*u + mat->b*u / ((u / (mat->E_0*eta*eta)+1.f) *
//        (u / (mat->E_0*eta*eta)+1.f)) *
//        (3.f*(u / (mat->E_0*eta*eta)+1.f) - 2.f) +
//        (mat->A + 2.f*mat->B*mu) / mat->rho_0  +  P_c / (rho*rho) *
//        (mat->a*rho + mat->b*rho / ((u / (mat->E_0*eta*eta)+1.f) *
//        (u / (mat->E_0*eta*eta)+1.f)));
//
//    c_c = max(c_c, mat->A / mat->rho_0);
//
//    // Expanded and hot
//    P_e = mat->a*density*u + (
//        mat->b * density * u / (u / (mat->E_0 * eta*eta) + 1.f)
//        + mat->A*mu * exp(-mat->beta * nu)
//        ) * exp(-mat->alpha * nu*nu);
//
//    c_e = (mat->a + mat->b / (u / (mat->E_0*eta*eta)+1.f) *
//        exp(-mat->beta*((1.f - eta)/eta)*((1.f - eta)/eta))
//        + 1.f)*P_e/rho + mat->A/mat->rho_0
//        *exp(-(mat->alpha*((1.f - eta)/eta)+mat->beta *
//        ((1.f - eta)/eta)*((1.f - eta)/eta)))*(1.f+mu/(eta*eta)
//        *(mat->alpha+2.f*mat->beta*((1.f - eta)/eta)-eta)) +
//        mat->b*rho*u/((u / (mat->E_0*eta*eta)+1.f)*
//        (u / (mat->E_0*eta*eta)+1.f)*eta*eta)*
//        exp(-mat->beta*((1.f - eta)/eta)*((1.f - eta)/eta))*
//        (2.f*mat->beta*((1.f - eta)/eta)*(u / (mat->E_0*eta*eta)+1.f) /
//         mat->rho_0 + 1.f/(mat->E_0*rho)*(2.f*u-P_e/rho));
//
//    // Condensed or cold state
//    if ((1.f < eta) || (u < mat->E_iv)) {
//        c = c_c;
//    }
//    // Expanded and hot state
//    else if ((eta < 1.f) && (mat->E_cv < u)) {
//        c = c_e;
//    }
//    // Hybrid state
//    else {
//		c = ((u - mat->E_iv)*c_e + (mat->E_cv - u)*c_c) /
//            (mat->E_cv - mat->E_iv);
//
//        c = max(c_c, mat->A / mat->rho0);
//    }
    float c;

    c = mat->c_TEMPORARY; /// VERY TEMPORARY!!!

    return c;
}

// gas_soundspeed_from_pressure
INLINE static float Til_soundspeed_from_pressure(
    float density, float P, struct Til_params *mat) {

    float c;

    c = mat->c_TEMPORARY; /// VERY TEMPORARY!!!

    return c;
}

#endif /* SWIFT_TILLOTSON_EQUATION_OF_STATE_H */


















