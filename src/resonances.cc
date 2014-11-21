/*
 *
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include "include/resonances.h"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_coupling.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <utility>
#include <vector>

#include "include/constants.h"
#include "include/decaymodes.h"
#include "include/distributions.h"
#include "include/fourvector.h"
#include "include/logging.h"
#include "include/macros.h"
#include "include/particledata.h"
#include "include/particles.h"
#include "include/particletype.h"
#include "include/processbranch.h"
#include "include/random.h"
#include "include/width.h"

namespace Smash {


double clebsch_gordan(const int j_a, const int j_b, const int j_c,
                      const int m_a, const int m_b, const int m_c) {
  const auto &log = logger<LogArea::Resonances>();
  double wigner_3j =  gsl_sf_coupling_3j(j_a, j_b, j_c, m_a, m_b, -m_c);
  double result = 0.;
  if (std::abs(wigner_3j) > really_small)
    result = std::pow(-1, (j_a-j_b+m_c)/2.) * std::sqrt(j_c + 1) * wigner_3j;

  log.debug("CG: ", result, " I1: ", j_a, " I2: ", j_b, " IR: ", j_c, " iz1: ", m_a,
            " iz2: ", m_b, " izR: ", m_c);

  return result;
}


/* two_to_one_formation -- only the resonance in the final state */
double two_to_one_formation(const ParticleType &type_particle_a,
                            const ParticleType &type_particle_b,
                            const ParticleType &type_resonance,
                            double mandelstam_s, double cm_momentum_squared) {
  /* Check for charge conservation */
  if (type_resonance.charge() != type_particle_a.charge()
                               + type_particle_b.charge()) {
    return 0.;
  }

  /* Check for baryon number conservation */
  if (type_resonance.baryon_number() != type_particle_a.baryon_number()
                                      + type_particle_b.baryon_number()) {
    return 0.;
  }

  /* Calculate partial width. */
  double srts = std::sqrt(mandelstam_s);
  float partial_width = type_resonance.get_partial_width(srts,
                        type_particle_a.pdgcode(), type_particle_b.pdgcode());
  if (partial_width <= 0.) {
    return 0.;
  }

  /* Calculate spin factor */
  const double spinfactor = (type_resonance.spin() + 1)
    / ((type_particle_a.spin() + 1) * (type_particle_b.spin() + 1));
  const int sym_factor = (type_particle_a.pdgcode()==type_particle_b.pdgcode())
                       ? 2 : 1;
  float resonance_width = type_resonance.total_width(srts);
  float resonance_mass = type_resonance.mass();
  /* Calculate resonance production cross section
   * using the Breit-Wigner distribution as probability amplitude.
   * See Eq. (176) in Buss et al., Physics Reports 512, 1 (2012). */
  return spinfactor * sym_factor * 4.0 * M_PI / cm_momentum_squared
         * breit_wigner(mandelstam_s, resonance_mass, resonance_width)
         * partial_width/resonance_width
         * hbarc * hbarc / fm2_mb;
}


/* Function for 1-dimensional GSL integration  */
void quadrature_1d(double (*integrand_function)(double, void*),
                     IntegrandParameters *parameters,
                     double lower_limit, double upper_limit,
                     double *integral_value, double *integral_error) {
  gsl_integration_workspace *workspace
    = gsl_integration_workspace_alloc(1000);
  gsl_function integrand;
  integrand.function = integrand_function;
  integrand.params = parameters;
  size_t subintervals_max = 100;
  int gauss_points = 2;
  double accuracy_absolute = 1.0e-6;
  double accuracy_relative = 1.0e-4;

  gsl_integration_qag(&integrand, lower_limit, upper_limit,
                      accuracy_absolute, accuracy_relative,
                      subintervals_max, gauss_points, workspace,
                      integral_value, integral_error);

  gsl_integration_workspace_free(workspace);
}

/* Spectral function of the resonance */
double spectral_function(double resonance_mass, double resonance_pole,
                         double resonance_width) {
  /* breit_wigner is essentially pi * mass * width * spectral function
   * (mass^2 is called mandelstam_s in breit_wigner)
   */
  return breit_wigner(resonance_mass * resonance_mass,
                      resonance_pole, resonance_width)
         / M_PI / resonance_mass / resonance_width;
}

/* Spectral function integrand for GSL integration */
double spectral_function_integrand(double resonance_mass,
                                   void *parameters) {
  IntegrandParameters *params
    = reinterpret_cast<IntegrandParameters*>(parameters);
  double resonance_pole_mass = params->type->mass();
  double stable_mass = params->m2;
  double mandelstam_s = params->s;
  double resonance_width = params->type->total_width(resonance_mass);

  /* center-of-mass momentum of final state particles */
  if (mandelstam_s - (stable_mass + resonance_mass)
      * (stable_mass + resonance_mass) > 0.0) {
    double cm_momentum_final
      = std::sqrt((mandelstam_s - (stable_mass + resonance_mass)
              * (stable_mass + resonance_mass))
             * (mandelstam_s - (stable_mass - resonance_mass)
                * (stable_mass - resonance_mass))
             / (4 * mandelstam_s));

    /* Integrand is the spectral function weighted by the
     * CM momentum of final state
     * In addition, dm^2 = 2*m*dm
     */
    return spectral_function(resonance_mass, resonance_pole_mass,
                             resonance_width)
           * cm_momentum_final
           * 2 * resonance_mass;
  } else {
    return 0.0;
  }
}

/* Resonance mass sampling for 2-particle final state */
double sample_resonance_mass(const ParticleType &type_resonance,
                             const ParticleType &type_stable,
                             const float cms_energy) {
  /* Define distribution parameters */
  float mass_stable = type_stable.mass();
  IntegrandParameters params = {&type_resonance, mass_stable,
                                cms_energy * cms_energy};

  /* sample resonance mass from the distribution
   * used for calculating the cross section
   */
  double mass_resonance = 0.0, random_number = 1.0;
  double distribution_max
    = spectral_function_integrand(params.type->mass(), &params);
  double distribution_value = 0.0;
  while (random_number > distribution_value) {
    random_number = Random::uniform(0.0, distribution_max);
    mass_resonance = Random::uniform(type_resonance.minimum_mass(),
                                     cms_energy - mass_stable);
    distribution_value = spectral_function_integrand(mass_resonance, &params);
  }

  return mass_resonance;
}


}  // namespace Smash
