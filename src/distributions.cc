/*
 *
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 *
 */
#include "include/distributions.h"

#include <stdio.h>
#include <stdint.h>
#include <gsl/gsl_sf_bessel.h>

#include "include/Box.h"
#include "include/constants.h"
#include "include/macros.h"
#include "include/ParticleType.h"

/* Breit-Wigner distribution for calculating resonance
 * production probability
 */
double breit_wigner(const double mandelstam_s, const float resonance_mass,
                    const float resonance_width) {
  return mandelstam_s * resonance_width * resonance_width
          / ((mandelstam_s - resonance_mass * resonance_mass)
             * (mandelstam_s - resonance_mass * resonance_mass)
             + mandelstam_s * resonance_width * resonance_width);
}

/* density_integrand - Maxwell-Boltzmann distribution */
double inline density_integrand(double momentum, double temp, double mass) {
  return 4 * M_PI * momentum * momentum
    * exp(- sqrt(momentum * momentum + mass * mass) / temp);
}

/* sample_momenta - return thermal momenta */
double sample_momenta(const Box &box, const ParticleType &type) {
  double momentum_radial, momentum_average, momentum_min, momentum_max;
  double probability = 0, probability_max, probability_random = 1;

  /* Maxwell-Boltzmann average E <E>=3T + m * K_1(m/T) / K_2(m/T) */
  momentum_average = sqrt((3 * box.temperature()
    + type.mass() * gsl_sf_bessel_K1(type.mass() / box.temperature())
                  / gsl_sf_bessel_Kn(2, type.mass() / box.temperature()))
    * (3 * box.temperature()
    + type.mass() * gsl_sf_bessel_K1(type.mass() / box.temperature())
                  / gsl_sf_bessel_Kn(2, type.mass() / box.temperature()))
    - type.mass() * type.mass());

  momentum_min = type.mass();
  momentum_max = 50.0 * box.temperature();
  /* double the massless peak value to be above maximum of the distribution */
  probability_max = 2 * density_integrand(momentum_average, box.temperature(),
    type.mass());

  /* sample by rejection method: (see numerical recipes for more efficient)
   * random momenta and random probability need to be below the distribution
   */
  while (probability_random > probability) {
    momentum_radial = (momentum_max - momentum_min) * drand48() + momentum_min;
    momentum_radial = sqrt(momentum_radial * momentum_radial
      - type.mass() * type.mass());
    probability = density_integrand(momentum_radial, box.temperature(),
      type.mass());
    probability_random = probability_max * drand48();
  }

  return momentum_radial;
}
