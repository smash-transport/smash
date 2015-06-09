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

#include "include/distributions.h"
#include "include/kinematics.h"

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

/* Spectral function of the resonance */
double spectral_function(double resonance_mass, double resonance_pole,
                         double resonance_width) {
  /* breit_wigner is essentially pi * mass * width * spectral function
   * (mass^2 is called mandelstam_s in breit_wigner)
   */
  return breit_wigner(resonance_mass * resonance_mass, resonance_pole,
                      resonance_width) /
         (M_PI * resonance_mass * resonance_width);
}

/* Integrand for spectral-function integration */
double spectral_function_integrand(double resonance_mass, double srts,
                                   double stable_mass, ParticleTypePtr type) {
  const double resonance_width = type->total_width(resonance_mass);

  if (srts < stable_mass + resonance_mass
      || resonance_width < really_small) {
    return 0.0;
  }

  /* Integrand is the spectral function weighted by the CM momentum of the
    * final state. In addition, dm^2 = 2*m*dm. */
  const double res_pole_mass = type->mass();
  return spectral_function(resonance_mass, res_pole_mass, resonance_width)
          * pCM(srts, stable_mass, resonance_mass)
          * 2 * resonance_mass;
}

/* Resonance mass sampling for 2-particle final state */
float sample_resonance_mass(const ParticleTypePtr type_resonance,
                            const ParticleTypePtr type_stable,
                            const double cms_energy) {
  /* Define distribution parameters */
  const float mass_stable = type_stable->mass();

  /* Sample resonance mass from the distribution
   * used for calculating the cross section. */
  float mass_resonance = 0.;
  float maximum_mass = std::nextafter(static_cast<float>(cms_energy -
                                                         mass_stable), 0.f);
  double random_number = 1.0;
  double distribution_max = spectral_function_integrand(type_resonance->mass(),
                                                        cms_energy, mass_stable,
                                                        type_resonance);
  double distribution_value = 0.0;
  while (random_number > distribution_value) {
    random_number = Random::uniform(0.0, distribution_max);
    mass_resonance = Random::uniform(type_resonance->minimum_mass(),
                                     maximum_mass);
    distribution_value = spectral_function_integrand(mass_resonance, cms_energy,
                                                     mass_stable,
                                                     type_resonance);
  }

  return mass_resonance;
}


}  // namespace Smash
