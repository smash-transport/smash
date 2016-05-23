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

/* relativistic Breit-Wigner distribution */
float breit_wigner(const float m, const float pole, const float width) {
  const float msqr = m*m;
  const float dmsqr = msqr - pole*pole;
  return 2.*msqr*width / (M_PI * (dmsqr*dmsqr + msqr*width*width));
}

/* non-relativistic Breit-Wigner distribution */
float breit_wigner_nonrel(float m, float pole, float width) {
  return cauchy(m, pole, width/2.);
}

float cauchy(float x, float pole, float width) {
  const float dm = x - pole;
  return width / (M_PI * (dm*dm + width*width));
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
  return 1.0/(std::exp((sqrt(momentum_radial*momentum_radial+mass*mass)
                        -baryon_chemical_potential)/temperature) + lam);
}


/* sample_momenta - return thermal momenta */
double sample_momenta(const double temperature, const double mass) {
  const auto &log = logger<LogArea::Distributions>();
  log.debug("Sample momenta with mass ", mass, " and T ", temperature);
  /* Maxwell-Boltzmann average E <E>=3T + m * K_1(m/T) / K_2(m/T) */
  float energy_average;
  if (mass > 0.) {
    // massive particles
    const float m_over_T = mass / temperature;
    energy_average = 3 * temperature
                     + mass * gsl_sf_bessel_K1(m_over_T)
                            / gsl_sf_bessel_Kn(2, m_over_T);
  } else {
    // massless particles
    energy_average = 3 * temperature;
  }
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
double sample_momenta_from_thermal(const double temperature,
                                   const double mass) {
  const auto &log = logger<LogArea::Distributions>();
  log.debug("Sample momenta with mass ", mass, " and T ", temperature);
  float momentum_radial, energy;
  // when temperature/mass
  if ( temperature > 0.6f*mass ) {
    while ( true ) {
      const float a = -std::log(Random::canonical_nonzero());
      const float b = -std::log(Random::canonical_nonzero());
      const float c = -std::log(Random::canonical_nonzero());
      momentum_radial = temperature * (a + b + c);
      energy = sqrt(momentum_radial * momentum_radial + mass * mass);
      if (Random::canonical() < exp((momentum_radial-energy)/temperature)) {
        break;
      }
    }
  } else {
    while ( true ) {
      const float r0 = Random::canonical();
      const float I1 = mass*mass;
      const float I2 = 2.0*mass*temperature;
      const float I3 = 2.0*temperature*temperature;
      const float Itot = I1 + I2 + I3;
      float K;
      if ( r0 < I1/Itot ) {
        const float r1 = Random::canonical_nonzero();
        K = -temperature*std::log(r1);
      } else if ( r0 < (I1+I2)/Itot ) {
        const float r1 = Random::canonical_nonzero();
        const float r2 = Random::canonical_nonzero();
        K = -temperature*std::log(r1*r2);
      } else {
        const float r1 = Random::canonical_nonzero();
        const float r2 = Random::canonical_nonzero();
        const float r3 = Random::canonical_nonzero();
        K = -temperature*std::log(r1*r2*r3);
      }
      energy = K + mass;
      momentum_radial = sqrt((energy + mass)*(energy - mass));
      if (Random::canonical() < momentum_radial/energy) {
        break;
      }
    }
  }
  return momentum_radial;
}

}  // namespace Smash
