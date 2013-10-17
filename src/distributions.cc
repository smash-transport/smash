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

#include <gsl/gsl_sf_bessel.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "include/outputroutines.h"

/* Breit-Wigner distribution for calculating resonance
 * production probability
 */
double breit_wigner(const double mandelstam_s, float resonance_mass,
                    float resonance_width) {
  return mandelstam_s * resonance_width * resonance_width
          / ((mandelstam_s - resonance_mass * resonance_mass)
             * (mandelstam_s - resonance_mass * resonance_mass)
             + mandelstam_s * resonance_width * resonance_width);
}

/* density_integrand - Maxwell-Boltzmann distribution */
double inline density_integrand(const double &momentum, const double &temp,
  const double &mass) {
  return 4 * M_PI * momentum * momentum
    * exp(- sqrt(momentum * momentum + mass * mass) / temp);
}

/* sample_momenta - return thermal momenta */
double sample_momenta(const double &temperature, const double &mass) {
  double momentum_radial, momentum_average, momentum_min, momentum_max;
  double probability = 0, probability_max, probability_random = 1;

  printd("Sample momenta with mass %g and T %g\n", mass, temperature);
  /* Maxwell-Boltzmann average E <E>=3T + m * K_1(m/T) / K_2(m/T) */
  momentum_average = sqrt((3 * temperature
    + mass * gsl_sf_bessel_K1(mass / temperature)
                  / gsl_sf_bessel_Kn(2, mass / temperature))
    * (3 * temperature + mass * gsl_sf_bessel_K1(mass / temperature)
                  / gsl_sf_bessel_Kn(2, mass / temperature)) - mass * mass);

  momentum_min = mass;
  momentum_max = 50.0 * temperature;
  /* double the massless peak value to be above maximum of the distribution */
  probability_max = 2 * density_integrand(momentum_average, temperature, mass);

  /* sample by rejection method: (see numerical recipes for more efficient)
   * random momenta and random probability need to be below the distribution
   */
  while (probability_random > probability) {
    momentum_radial = (momentum_max - momentum_min) * drand48() + momentum_min;
    momentum_radial = sqrt(momentum_radial * momentum_radial - mass * mass);
    probability = density_integrand(momentum_radial, temperature, mass);
    probability_random = probability_max * drand48();
  }

  return momentum_radial;
}
