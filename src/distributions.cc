/*
 *
 *    Copyright (c) 2013-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include "include/distributions.h"

#include <gsl/gsl_sf_bessel.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "include/constants.h"
#include "include/logging.h"
#include "include/random.h"

namespace Smash {

/* Breit-Wigner distribution for calculating resonance production probability */
float breit_wigner(const double mandelstam_s, const float resonance_mass,
    const float resonance_width) {
    const double A = mandelstam_s * resonance_width * resonance_width;
    const double B = (mandelstam_s - resonance_mass * resonance_mass);
    return A / (B*B + A);
}

/** woods-saxon distribution function  */
double woods_saxon_dist_func(const double r,  const double radius,
        const double diffusion) {
    return 1.0/(std::exp((r-radius)/diffusion)+1.0);
}

/* density_integrand - Maxwell-Boltzmann distribution */
double density_integrand(const double energy, const double momentum_sqr,
                         const double temperature) {
  return 4.0 * M_PI * momentum_sqr * exp(-energy / temperature);
}

/* General Juttner distribution
 * lam = 1: to Fermion-Dirac distribution
 * lam = 0: to Thermal distribution
 * lam =-1: to Bose-Einstein distribution */
double juttner_distribution_func(const double momentum_radial,
        const double mass, const double temperature, const double
        baryon_chemical_potential, const double lam) {
    return 1.0/(std::exp((sqrt(momentum_radial*
                     momentum_radial+mass*mass)-baryon_chemical_potential)/
                temperature) + lam);
}


/* sample_momenta - return thermal momenta */
double sample_momenta(const double temperature, const double mass) {
  const auto &log = logger<LogArea::Distributions>();
  log.debug("Sample momenta with mass ", mass, " and T ", temperature);
  /* Maxwell-Boltzmann average E <E>=3T + m * K_1(m/T) / K_2(m/T) */
  const float m_over_T = mass / temperature;
  const float energy_average = 3 * temperature
                             + mass * gsl_sf_bessel_K1(m_over_T)
                                    / gsl_sf_bessel_Kn(2, m_over_T);
  const float momentum_average_sqr = (energy_average - mass) *
                                     (energy_average + mass);
  const float energy_min = mass;
  const float energy_max = 50.0f * temperature;
  /* double the massless peak value to be above maximum of the distribution */
  const float probability_max = 2.0f * density_integrand(energy_average,
                                                         momentum_average_sqr,
                                                         temperature);

  /* sample by rejection method: (see numerical recipes for more efficient)
   * random momenta and random probability need to be below the distribution */
  float momentum_radial_sqr, probability;
  do {
    float energy = Random::uniform(energy_min, energy_max);
    momentum_radial_sqr = (energy - mass) * (energy + mass);
    probability = density_integrand(energy, momentum_radial_sqr, temperature);
  } while (Random::uniform(0.f, probability_max) > probability);

  return std::sqrt(momentum_radial_sqr);
}

/* A much faster sampler for thermal distribution from Scott
* (see \iref{Pratt:2014vja}) APPENDIX: ALGORITHM FOR GENERATING PARTICLES
* math trick: for \f$ x^{n-1}e^{-x} \f$ distribution, sample x by:
* \f$ x = -ln(r_1 r_2 r_3 ... r_n) \f$
* where \f$ r_i \f$ are uniform random numbers between [0,1)
* for \f$ T/m > 0.6 \f$: \f$ p^2 e^{-E/T} = p^2 e^{-p/T} * e^{(p-E)/T} \f$,
* where \f$ e^{(p-E)/T}\f$ is used as rejection weight.
* Since \f$T/m > 0.6 \f$, \f$ e^{(p-E)/T}\f$ is close to 1.
* for \f$ T/m < 0.6 \f$, there are many rejections
* another manipulation is used:
* \f$ p^2 e^{-E/T} dp = dE \frac{E}{p} p^2 e^{-E/T} \f$
* \f$ = dK \frac{p}{E} (K+m)^2 e^{-K/T} e^{-m/T} \f$
* \f$ = dK (K^2 + 2mK + m^2) e^{-K/T} \frac{p}{E}\f$
*  where \frac{p}{E} is used as rejection weight.
* return: themal momenta */

double sample_momenta_from_thermal(const double temperature, const double mass) {
    const auto &log = logger<LogArea::Distributions>();
    log.debug("Sample momenta with mass ", mass, " and T ", temperature);
    float momentum_radial, energy;
    float r0, r1, r2, r3, a, b, c;
    float K, I1, I2, I3, Itot;
    // when temperature/mass
    if ( temperature/mass > 0.6f ) {
        while ( true ) {
            r1 = Random::canonical();
            r2 = Random::canonical();
            r3 = Random::canonical();
            a = -std::log(r1);
            b = -std::log(r2);
            c = -std::log(r3);
            momentum_radial = temperature * (a + b + c);
            energy = sqrt(momentum_radial * momentum_radial + mass * mass);
            if ( Random::canonical() <
                    exp((momentum_radial-energy)/temperature) ) {
                break;
            }
        }
    } else {
        while ( true ) {
            r0 = Random::canonical();
            I1 = mass*mass;
            I2 = 2.0*mass*temperature;
            I3 = 2.0*temperature*temperature;
            Itot = I1 + I2 + I3;
            if ( r0 < I1/Itot ) {
                r1 = Random::canonical();
                K = -temperature*std::log(r1);
            } else if ( r0 < (I1+I2)/Itot ) {
                r1 = Random::canonical();
                r2 = Random::canonical();
                K = -temperature*std::log(r1*r2);
            } else {
                r1 = Random::canonical();
                r2 = Random::canonical();
                r3 = Random::canonical();
                K = -temperature*std::log(r1*r2*r3);
            }
            energy = K + mass;
            momentum_radial = sqrt(energy*energy - mass*mass);
            r0 = Random::canonical();
            if ( r0 < momentum_radial/energy ) break;
        }
    }
    return momentum_radial;
}
}  // namespace Smash
