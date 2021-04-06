/**
 * @file blood_consts.h
 * @author Ilya Kapitonau <ilya.th.kapitonov@gmail.com>
 * @brief Provides blood constants for Carreau fluid model.
 *
 * According to Carreau model, fluid viscosity is calculated by the following
 * equation: mu = mu_inf + (mu_0 - mu_inf) * (1 + (b * shear_rate)^2)^((n - 1) / 2).
 */

#ifndef BLOOD_CONSTS_H_
#define BLOOD_CONSTS_H_

/**
 * @brief Viscosity at zero shear rate (Pa.s).
 */
extern const double MU_0;

/**
 * @brief Viscosity at infinite shear rate (Pa.s).
 */
extern const double MU_INF;

/**
 * @brief Relaxation time (s).
 */
extern const double B;

/**
 * @brief Power index.
 */
extern const double N;

#endif	// BLOOD_CONSTS_H_
