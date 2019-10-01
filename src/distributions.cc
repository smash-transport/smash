/*
 *
 *    Copyright (c) 2013-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include "smash/distributions.h"

#include <gsl/gsl_sf_bessel.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "smash/constants.h"
#include "smash/logging.h"
#include "smash/random.h"

namespace smash {
inline constexpr int Distributions = LogArea::Distributions::id;

// Relativistic Breit-Wigner distribution
double breit_wigner(const double m, const double pole, const double width) {
  const double msqr = m * m;
  const double dmsqr = msqr - pole * pole;
  return 2. * msqr * width / (M_PI * (dmsqr * dmsqr + msqr * width * width));
}

// Non-relativistic Breit-Wigner distribution
double breit_wigner_nonrel(double m, double pole, double width) {
  return cauchy(m, pole, width / 2.);
}

// Cauchy distribution used in non-relativistic Breit-Wigner above
double cauchy(double x, double pole, double width) {
  const double dm = x - pole;
  return width / (M_PI * (dm * dm + width * width));
}

// density_integrand - Maxwell-Boltzmann distribution
double density_integrand(const double energy, const double momentum_sqr,
                         const double temperature) {
  return 4.0 * M_PI * momentum_sqr * std::exp(-energy / temperature);
}

double density_integrand_mass(const double energy, const double momentum_sqr,
                              const double temperature) {
  return momentum_sqr * std::sqrt(momentum_sqr) *
         std::exp(-energy / temperature);
}

double density_integrand_1M_IC(const double energy, const double momentum_sqr,
                               const double temperature) {
  return ((3.0 / 20.0) * (momentum_sqr / (temperature * temperature)) -
          (6.0 / 5.0) * (energy / temperature) + (14.0 / 5.0)) *
         std::exp(-energy / temperature) * momentum_sqr;
}

double density_integrand_2M_IC(const double energy, const double momentum_sqr,
                               const double temperature) {
  return (0.75 + 0.125 * (momentum_sqr / (temperature * temperature)) -
          (1.0 / 30.0) * (momentum_sqr * energy /
                          (temperature * temperature * temperature)) +
          (1.0 / 480.0) *
              (momentum_sqr * momentum_sqr /
               (temperature * temperature * temperature * temperature))) *
         std::exp(-energy / temperature) * momentum_sqr;
}

double juttner_distribution_func(const double momentum_radial,
                                 const double mass, const double temperature,
                                 const double baryon_chemical_potential,
                                 const double lam) {
  return 1.0 /
         (std::exp((std::sqrt(momentum_radial * momentum_radial + mass * mass) -
                    baryon_chemical_potential) /
                   temperature) +
          lam);
}

double sample_momenta_non_eq_mass(const double temperature, const double mass) {
  logg[Distributions].debug("Sample momenta with mass ", mass, " and T ",
                            temperature);

  /* Calculate range on momentum values to use
   * ideally 0.0 and as large as possible but we want to be efficient! */
  const double mom_min = 0.0;
  const double mom_max =
      std::sqrt(50. * 50. * temperature * temperature - mass * mass);
  /* Calculate momentum and energy values that will give
   * maxima of density_integrand_mass, verified by differentiation */
  const double p_non_eq_sq =
      0.5 * (9 * temperature * temperature +
             temperature *
                 std::sqrt(81 * temperature * temperature + 36 * mass * mass));
  const double e_non_eq = std::sqrt(p_non_eq_sq + mass * mass);
  const double probability_max =
      2. * density_integrand_mass(e_non_eq, p_non_eq_sq, temperature);

  /* sample by rejection method: (see numerical recipes for more efficient)
   * random momenta and random probability need to be below the distribution */
  double energy, momentum_radial, probability;
  do {
    // sample uniformly in momentum, DONT sample uniformly in energy!
    momentum_radial = random::uniform(mom_min, mom_max);
    // Energy by on-shell condition
    energy = std::sqrt(momentum_radial * momentum_radial + mass * mass);
    probability = density_integrand_mass(
        energy, momentum_radial * momentum_radial, temperature);
  } while (random::uniform(0., probability_max) > probability);

  return momentum_radial;
}

double sample_momenta_1M_IC(const double temperature, const double mass) {
  logg[Distributions].debug("Sample momenta with mass ", mass, " and T ",
                            temperature);
  // Maxwell-Boltzmann average E <E>=3T + m * K_1(m/T) / K_2(m/T)
  double energy_average;
  if (mass > 0.) {
    // massive particles
    const double m_over_T = mass / temperature;
    energy_average = 3 * temperature + mass * gsl_sf_bessel_K1(m_over_T) /
                                           gsl_sf_bessel_Kn(2, m_over_T);
  } else {
    // massless particles
    energy_average = 3 * temperature;
  }
  const double momentum_average_sqr =
      (energy_average - mass) * (energy_average + mass);
  const double energy_min = mass;
  const double energy_max = 50. * temperature;
  /* 16 * the massless peak value to be well above maximum of the
   * distribution */
  const double probability_max =
      16. * density_integrand_1M_IC(energy_average, momentum_average_sqr,
                                    temperature);

  /* sample by rejection method: (see numerical recipes for more efficient)
   * random momenta and random probability need to be below the distribution */
  double momentum_radial_sqr, probability;
  do {
    double energy = random::uniform(energy_min, energy_max);
    momentum_radial_sqr = (energy - mass) * (energy + mass);
    probability =
        density_integrand_1M_IC(energy, momentum_radial_sqr, temperature);
  } while (random::uniform(0., probability_max) > probability);

  return std::sqrt(momentum_radial_sqr);
}

double sample_momenta_2M_IC(const double temperature, const double mass) {
  logg[Distributions].debug("Sample momenta with mass ", mass, " and T ",
                            temperature);
  /* Maxwell-Boltzmann average E <E>=3T + m * K_1(m/T) / K_2(m/T) */
  double energy_average;
  if (mass > 0.) {
    // massive particles
    const double m_over_T = mass / temperature;
    energy_average = 3 * temperature + mass * gsl_sf_bessel_K1(m_over_T) /
                                           gsl_sf_bessel_Kn(2, m_over_T);
  } else {
    // massless particles
    energy_average = 3 * temperature;
  }
  const double momentum_average_sqr =
      (energy_average - mass) * (energy_average + mass);
  const double energy_min = mass;
  const double energy_max = 50. * temperature;
  /* 16 * the massless peak value to be well above maximum of the
   * distribution */
  const double probability_max =
      16. * density_integrand_2M_IC(energy_average, momentum_average_sqr,
                                    temperature);

  /* sample by rejection method: (see numerical recipes for more efficient)
   * random momenta and random probability need to be below the distribution */
  double momentum_radial_sqr, probability;
  do {
    double energy = random::uniform(energy_min, energy_max);
    momentum_radial_sqr = (energy - mass) * (energy + mass);
    probability =
        density_integrand_2M_IC(energy, momentum_radial_sqr, temperature);
  } while (random::uniform(0., probability_max) > probability);

  return std::sqrt(momentum_radial_sqr);
}

double sample_momenta_from_thermal(const double temperature,
                                   const double mass) {
  logg[Distributions].debug("Sample momenta with mass ", mass, " and T ",
                            temperature);
  double momentum_radial, energy;
  // when temperature/mass
  if (temperature > 0.6 * mass) {
    while (true) {
      const double a = -std::log(random::canonical_nonzero());
      const double b = -std::log(random::canonical_nonzero());
      const double c = -std::log(random::canonical_nonzero());
      momentum_radial = temperature * (a + b + c);
      energy = std::sqrt(momentum_radial * momentum_radial + mass * mass);
      if (random::canonical() <
          std::exp((momentum_radial - energy) / temperature)) {
        break;
      }
    }
  } else {
    while (true) {
      const double r0 = random::canonical();
      const double I1 = mass * mass;
      const double I2 = 2.0 * mass * temperature;
      const double I3 = 2.0 * temperature * temperature;
      const double Itot = I1 + I2 + I3;
      double K;
      if (r0 < I1 / Itot) {
        const double r1 = random::canonical_nonzero();
        K = -temperature * std::log(r1);
      } else if (r0 < (I1 + I2) / Itot) {
        const double r1 = random::canonical_nonzero();
        const double r2 = random::canonical_nonzero();
        K = -temperature * std::log(r1 * r2);
      } else {
        const double r1 = random::canonical_nonzero();
        const double r2 = random::canonical_nonzero();
        const double r3 = random::canonical_nonzero();
        K = -temperature * std::log(r1 * r2 * r3);
      }
      energy = K + mass;
      momentum_radial = std::sqrt((energy + mass) * (energy - mass));
      if (random::canonical() < momentum_radial / energy) {
        break;
      }
    }
  }
  return momentum_radial;
}

double sample_momenta_IC_ES(const double temperature) {
  double momentum_radial;
  const double a = -std::log(random::canonical_nonzero());
  const double b = -std::log(random::canonical_nonzero());
  const double c = -std::log(random::canonical_nonzero());
  const double d = -std::log(random::canonical_nonzero());
  momentum_radial = (3.0 / 4.0) * temperature * (a + b + c + d);

  return momentum_radial;
}

}  // namespace smash
