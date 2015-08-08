/*
 *
 *    Copyright (c) 2013-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include "include/resonances.h"

#include <gsl/gsl_sf_coupling.h>

#include "include/kinematics.h"
#include "include/logging.h"
#include "include/random.h"

namespace Smash {


double clebsch_gordan(const int j_a, const int j_b, const int j_c,
                      const int m_a, const int m_b, const int m_c) {
  const auto &log = logger<LogArea::Resonances>();
  double wigner_3j =  gsl_sf_coupling_3j(j_a, j_b, j_c, m_a, m_b, -m_c);
  double result = 0.;
  if (std::abs(wigner_3j) > really_small)
    result = std::pow(-1, (j_a-j_b+m_c)/2.) * std::sqrt(j_c + 1) * wigner_3j;

  log.debug("CG: ", result, " I1: ", j_a, " I2: ", j_b, " IR: ", j_c,
            " iz1: ", m_a, " iz2: ", m_b, " izR: ", m_c);

  return result;
}

/* Integrand for spectral-function integration */
float spectral_function_integrand(float resonance_mass, float srts,
                                  float stable_mass, const ParticleType &type) {
  if (srts <= stable_mass + resonance_mass) {
    return 0.;
  }

  /* Integrand is the spectral function weighted by the CM momentum of the
   * final state. */
  return type.spectral_function(resonance_mass)
         * pCM(srts, stable_mass, resonance_mass);
}

/* Resonance mass sampling for 2-particle final state */
float sample_resonance_mass(const ParticleType &type_resonance,
                            const float mass_stable, const float cms_energy) {
  /* Sample resonance mass from the distribution
   * used for calculating the cross section. */
  float mass_resonance = 0.;
  float maximum_mass = std::nextafter(static_cast<float>(cms_energy -
                                                         mass_stable), 0.f);
  float random_number = 1.;
  float distribution_max = spectral_function_integrand(
                                  std::min(maximum_mass, type_resonance.mass()),
                                                        cms_energy, mass_stable,
                                                        type_resonance);
  float distribution_value = 0.;
  while (random_number > distribution_value) {
    random_number = Random::uniform(0.f, distribution_max);
    mass_resonance = Random::uniform(type_resonance.minimum_mass(),
                                     maximum_mass);
    distribution_value = spectral_function_integrand(mass_resonance, cms_energy,
                                                     mass_stable,
                                                     type_resonance);
  }

  if (distribution_value > distribution_max) {
    const auto &log = logger<LogArea::Resonances>();
    log.error("distribution value too large in sample_resonance_mass: ",
              distribution_value, " ", distribution_max, " ",
              type_resonance.pdgcode(), " ", mass_stable, " ", cms_energy, " ",
              mass_resonance);
  }

  return mass_resonance;
}


}  // namespace Smash
