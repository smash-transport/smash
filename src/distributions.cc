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
double breit_wigner(const double m, const double pole, const double width) {
  const double msqr = m*m;
  const double dmsqr = msqr - pole*pole;
  return 2.*msqr*width / (M_PI * (dmsqr*dmsqr + msqr*width*width));
}

/* non-relativistic Breit-Wigner distribution */
double breit_wigner_nonrel(double m, double pole, double width) {
  return cauchy(m, pole, width/2.);
}

double cauchy(double x, double pole, double width) {
  const double dm = x - pole;
  return width / (M_PI * (dm*dm + width*width));
}

/** woods-saxon distribution function  */
double woods_saxon_dist_func(const double r,  const double radius,
                             const double diffusion) {
  return 1.0/(std::exp((r-radius)/diffusion)+1.0);
}

/** density_integrand - Maxwell-Boltzmann distribution */
double density_integrand(const double energy, const double momentum_sqr,
                         const double temperature) {
  return 4.0 * M_PI * momentum_sqr * exp(-energy / temperature);
}

/** density_integrand - off_equilibrium distribution for massive particles
 * \f[f=pe^{-\frac{\sqrt{m^2+p^2}}{T_0}}\f] */
double density_integrand_mass(const double energy, const double momentum_sqr,
                         const double temperature) {
  return momentum_sqr * std::sqrt(momentum_sqr)* exp(-energy / temperature);
}

/** density integrand - 1M_IC massless particles for expanding metric initialization,
 * see \iref{Bazow:2016oky}*/
double density_integrand_1M_IC(const double energy, const double momentum_sqr,
                        const double temperature) {
  return ((3.0/20.0) * (momentum_sqr/(temperature*temperature)) - (6.0/5.0) *
    (energy/temperature) + (14.0/5.0))*exp(-energy / temperature)*momentum_sqr;
}

/** density integrand - 2M_IC massless particles for expanding metric initialization,
 * see \iref{Bazow:2016oky}*/
double density_integrand_2M_IC(const double energy, const double momentum_sqr,
                        const double temperature) {
  return (0.75 + 0.125 * (momentum_sqr / (temperature * temperature)) -
         (1.0/30.0) * (momentum_sqr * energy /
         (temperature * temperature * temperature)) +
         (1.0/480.0) * (momentum_sqr * momentum_sqr/
         (temperature * temperature * temperature * temperature)))
         *exp(-energy/temperature) * momentum_sqr;
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


/** sample_momenta via rejection method, according to non-equilibrium 
 * distribution \f[f=pe^{-\frac{\sqrt{m^2+p^2}}{T_0}}\f]
 */
double sample_momenta_non_eq_mass(const double temperature, const double mass) {
  const auto &log = logger<LogArea::Distributions>();
  log.debug("Sample momenta with mass ", mass, " and T ", temperature);

  // calculate range on momentum values to use
  // ideally 0.0 and as large as possible but we want to be efficient!
  const double mom_min = 0.0;
  const double mom_max = std::sqrt(50. * 50. * temperature * temperature
                                  - mass * mass);
  // calculate momentum and energy values that will give
  // maxima of density_integrand_mass, verified by differentiation
  const double p_non_eq_sq = 0.5*(9*temperature*temperature +
                  temperature*std::sqrt(81*temperature*temperature
                                        + 36*mass*mass));
  const double e_non_eq = std::sqrt(p_non_eq_sq + mass*mass);
  const double probability_max = 2.*
                  density_integrand_mass(e_non_eq, p_non_eq_sq, temperature);

  /* sample by rejection method: (see numerical recipes for more efficient)
   * random momenta and random probability need to be below the distribution */
  double energy, momentum_radial, probability;
  do {
    // sample uniformly in momentum, DONT sample uniformly in energy!
    momentum_radial = Random::uniform(mom_min, mom_max);
    // Energy by on-shell condition
    energy = std::sqrt(momentum_radial*momentum_radial + mass*mass);
    probability = density_integrand_mass(energy,
                          momentum_radial*momentum_radial, temperature);
  } while (Random::uniform(0., probability_max) > probability);

  return momentum_radial;
}

/** sample_momenta for the 1M IC condition - return thermal momenta */
double sample_momenta_IC_1M(const double temperature, const double mass) {
  const auto &log = logger<LogArea::Distributions>();
  log.debug("Sample momenta with mass ", mass, " and T ", temperature);
  /* Maxwell-Boltzmann average E <E>=3T + m * K_1(m/T) / K_2(m/T) */
  double energy_average;
  if (mass > 0.) {
    // massive particles
    const double m_over_T = mass / temperature;
    energy_average = 3 * temperature
                     + mass * gsl_sf_bessel_K1(m_over_T)
                            / gsl_sf_bessel_Kn(2, m_over_T);
  } else {
    // massless particles
    energy_average = 3 * temperature;
  }
  const double momentum_average_sqr = (energy_average - mass) *
                                     (energy_average + mass);
  const double energy_min = mass;
  const double energy_max = 50. * temperature;
  /* 16 * the massless peak value to be well above maximum of the
   * distribution */
  const double probability_max = 16. * density_integrand_1M_IC(energy_average,
                                                         momentum_average_sqr,
                                                         temperature);

  /* sample by rejection method: (see numerical recipes for more efficient)
   * random momenta and random probability need to be below the distribution */
  double momentum_radial_sqr, probability;
  do {
    double energy = Random::uniform(energy_min, energy_max);
    momentum_radial_sqr = (energy - mass) * (energy + mass);
    probability = density_integrand_1M_IC(energy,
                          momentum_radial_sqr, temperature);
  } while (Random::uniform(0., probability_max) > probability);

  return std::sqrt(momentum_radial_sqr);
}

/** sample_momenta for the 2M IC condition - return thermal momenta */
double sample_momenta_IC_2M(const double temperature, const double mass) {
  const auto &log = logger<LogArea::Distributions>();
  log.debug("Sample momenta with mass ", mass, " and T ", temperature);
  /* Maxwell-Boltzmann average E <E>=3T + m * K_1(m/T) / K_2(m/T) */
  double energy_average;
  if (mass > 0.) {
    // massive particles
    const double m_over_T = mass / temperature;
    energy_average = 3 * temperature
                     + mass * gsl_sf_bessel_K1(m_over_T)
                            / gsl_sf_bessel_Kn(2, m_over_T);
  } else {
    // massless particles
    energy_average = 3 * temperature;
  }
  const double momentum_average_sqr = (energy_average - mass) *
                                     (energy_average + mass);
  const double energy_min = mass;
  const double energy_max = 50. * temperature;
  /* 16 * the massless peak value to be well above maximum of the
   * distribution */
  const double probability_max = 16. * density_integrand_2M_IC(energy_average,
                                                         momentum_average_sqr,
                                                         temperature);

  /* sample by rejection method: (see numerical recipes for more efficient)
   * random momenta and random probability need to be below the distribution */
  double momentum_radial_sqr, probability;
  do {
    double energy = Random::uniform(energy_min, energy_max);
    momentum_radial_sqr = (energy - mass) * (energy + mass);
    probability = density_integrand_2M_IC(energy, momentum_radial_sqr,
                                          temperature);
  } while (Random::uniform(0., probability_max) > probability);

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
  double momentum_radial, energy;
  // when temperature/mass
  if ( temperature > 0.6*mass ) {
    while ( true ) {
      const double a = -std::log(Random::canonical_nonzero());
      const double b = -std::log(Random::canonical_nonzero());
      const double c = -std::log(Random::canonical_nonzero());
      momentum_radial = temperature * (a + b + c);
      energy = sqrt(momentum_radial * momentum_radial + mass * mass);
      if (Random::canonical() < exp((momentum_radial-energy)/temperature)) {
        break;
      }
    }
  } else {
    while ( true ) {
      const double r0 = Random::canonical();
      const double I1 = mass*mass;
      const double I2 = 2.0*mass*temperature;
      const double I3 = 2.0*temperature*temperature;
      const double Itot = I1 + I2 + I3;
      double K;
      if ( r0 < I1/Itot ) {
        const double r1 = Random::canonical_nonzero();
        K = -temperature*std::log(r1);
      } else if ( r0 < (I1+I2)/Itot ) {
        const double r1 = Random::canonical_nonzero();
        const double r2 = Random::canonical_nonzero();
        K = -temperature*std::log(r1*r2);
      } else {
        const double r1 = Random::canonical_nonzero();
        const double r2 = Random::canonical_nonzero();
        const double r3 = Random::canonical_nonzero();
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
/// sample momenta according to the momentum distribution in ref{Bazow:2016oky}
double sample_momenta_IC_ES(const double temperature) {
  double momentum_radial;
  const double a = -std::log(Random::canonical_nonzero());
  const double b = -std::log(Random::canonical_nonzero());
  const double c = -std::log(Random::canonical_nonzero());
  const double d = -std::log(Random::canonical_nonzero());
  momentum_radial = (3.0/4.0)*temperature * (a + b + c + d);

  return momentum_radial;
}

}  // namespace Smash
