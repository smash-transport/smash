/*
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_DISTRIBUTIONS_H_
#define SRC_INCLUDE_DISTRIBUTIONS_H_

#include <gsl/gsl_sf_bessel.h>
#include <cmath>

#include "include/constants.h"

namespace Smash {

/* Breit-Wigner distribution for calculating resonance
 * production probability
 */
double breit_wigner(double mandelstam_s, float resonance_mass,
                          float resonance_width);

/* density_integrand - Maxwell-Boltzmann distribution */
double inline density_integrand(const double &momentum, const double &temp,
  const double &mass);

/* sample_momenta - return thermal momenta */
double sample_momenta(const double &temp, const double &mass);

/* return number density for given mass and temperature */
inline double number_density_maxwellboltzmann(double mass, double temperature) {
  /*
   * The particle number depends on distribution function
   * (assumes Maxwell-Boltzmann):
   * Volume m^2 T BesselK[2, m/T] / (2\pi^2)
   */
  return mass * mass * temperature * gsl_sf_bessel_Knu(2, mass / temperature)
    * 0.5 * M_1_PI * M_1_PI / hbarc / hbarc / hbarc;
}

}  // namespace Smash

#endif  // SRC_INCLUDE_DISTRIBUTIONS_H_
