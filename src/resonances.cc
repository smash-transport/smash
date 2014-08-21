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
#include "include/macros.h"
#include "include/outputroutines.h"
#include "include/particles.h"
#include "include/particledata.h"
#include "include/particletype.h"
#include "include/processbranch.h"
#include "include/random.h"
#include "include/width.h"

namespace Smash {

/* Parameters for spectral function integration via GSL. */
struct IntegrandParameters {
  const ParticleType *type;
  double m2;
  double s;
};


double clebsch_gordan(const int j1, const int j2, const int j3,
                      const int m1, const int m2, const int m3) {
  double wigner_3j =  gsl_sf_coupling_3j(j1, j2, j3, m1, m2, -m3);
  double result = 0.;
  if (std::abs(wigner_3j) > really_small)
    result = std::pow(-1, (j1-j2+m3)/2.) * std::sqrt(j3 + 1) * wigner_3j;

  printd("CG: %g I1: %i I2: %i IR: %i iz1: %i iz2: %i izR: %i \n",
         result, j1, j2, j3, m1, m2, m3);

  return result;
}


/* Calculate isospin Clebsch-Gordan coefficient for two particles p1 and p2
 * coupling to a total isospin (I_tot, I_z). */
inline double isospin_clebsch_gordan(const ParticleType p1,
                                     const ParticleType p2,
                                     const int I_tot, const int I_z) {
  return clebsch_gordan (p1.isospin(), p2.isospin(), I_tot,
                         p1.isospin3(), p2.isospin3(), I_z);
}

/* Calculate isospin Clebsch-Gordan coefficient for two particles p1 and p2
 * coupling to a resonance Res. */
inline double isospin_clebsch_gordan(const ParticleType p1,
                                     const ParticleType p2,
                                     const ParticleType Res) {
  return clebsch_gordan (p1.isospin(), p2.isospin(), Res.isospin(),
                         p1.isospin3(), p2.isospin3(), Res.isospin3());
}

/**
 * Function for 1-dimensional GSL integration.
 *
 * \param[in] integrand_function Function of 1 variable to be integrated over.
 * \param[in] parameters Container for possible parameters
 * needed by the integrand.
 * \param[in] lower_limit Lower limit of the integral.
 * \param[in] upper_limit Upper limit of the integral.
 * \param[out] integral_value Result of integration.
 * \param[out] integral_error Uncertainty of the result.
 */
static void quadrature_1d(double (*integrand_function)(double, void *),
                          IntegrandParameters *parameters, double lower_limit,
                          double upper_limit, double *integral_value,
                          double *integral_error);


/* two_to_one_formation -- only the resonance in the final state */
double two_to_one_formation(const ParticleType &type_particle1,
                            const ParticleType &type_particle2,
                            const ParticleType &type_resonance,
                            double mandelstam_s, double cm_momentum_squared) {
  /* Check for charge conservation */
  if (type_resonance.charge() != type_particle1.charge()
                                 + type_particle2.charge())
    return 0.0;

  /* Check for baryon number conservation */
  if (type_particle1.spin() % 2 != 0 || type_particle2.spin() % 2 != 0) {
    /* Step 1: We must have fermion */
    if (type_resonance.spin() % 2 == 0) {
      return 0.0;
    }
    /* Step 2: We must have antiparticle for antibaryon
     * (and non-antiparticle for baryon)
     */
    if (type_particle1.pdgcode().baryon_number() != 0
        && (type_particle1.pdgcode().baryon_number()
            != type_resonance.pdgcode().baryon_number())) {
      return 0.0;
    } else if (type_particle2.pdgcode().baryon_number() != 0
        && (type_particle2.pdgcode().baryon_number()
        != type_resonance.pdgcode().baryon_number())) {
      return 0.0;
    }
  }

  double isospin_factor
    = isospin_clebsch_gordan(type_particle1, type_particle2, type_resonance);

  /* If Clebsch-Gordan coefficient is zero, don't bother with the rest */
  if (std::abs(isospin_factor) < really_small)
    return 0.0;

  /* Check the decay modes of this resonance */
  const std::vector<DecayBranch> decaymodes
    = DecayModes::find(type_resonance.pdgcode()).decay_mode_list();
  bool not_enough_energy = false;
  /* Detailed balance required: Formation only possible if
   * the resonance can decay back to these particles
   */
  bool not_balanced = true;
  for (const auto &mode : decaymodes) {
    size_t decay_particles = mode.pdg_list().size();
    if ( decay_particles > 3 ) {
      printf("Warning: Not a 1->2 or 1->3 process!\n");
      printf("Number of decay particles: %zu \n", decay_particles);
    } else {
      /* There must be enough energy to produce all decay products */
      if (std::sqrt(mandelstam_s) < mode.threshold())
        not_enough_energy = true;
      /* Initial state is also a possible final state;
       * weigh the cross section with the ratio of this branch
       * XXX: For now, assuming only 2-particle initial states
       */
      if (decay_particles == 2
          && ((mode.pdg_list().at(0) == type_particle1.pdgcode()
               && mode.pdg_list().at(1) == type_particle2.pdgcode())
              || (mode.pdg_list().at(0) == type_particle2.pdgcode()
                  && mode.pdg_list().at(1) == type_particle1.pdgcode()))
          && (mode.weight() > 0.0))
        not_balanced = false;
    }
  }
  if (not_enough_energy || not_balanced) {
    return 0.0;
  }

  /* Calculate spin factor */
  const double spinfactor = (type_resonance.spin() + 1)
    / ((type_particle1.spin() + 1) * (type_particle2.spin() + 1));
  float resonance_width = type_resonance.total_width(std::sqrt(mandelstam_s));
  float resonance_mass = type_resonance.mass();
  /* Calculate resonance production cross section
   * using the Breit-Wigner distribution as probability amplitude
   * See Eq. (176) in Buss et al., Physics Reports 512, 1 (2012)
   */
  return isospin_factor * isospin_factor * spinfactor
         * 4.0 * M_PI / cm_momentum_squared
         * breit_wigner(mandelstam_s, resonance_mass, resonance_width)
         * hbarc * hbarc / fm2_mb;
}

/* two_to_two_formation -- resonance and another particle in final state */
size_t two_to_two_formation(const ParticleType &type_particle1,
                            const ParticleType &type_particle2,
                            const ParticleType &type_resonance,
                            double mandelstam_s, double cm_momentum_squared,
                            std::vector<ProcessBranch> *process_list) {
  size_t number_of_processes = 0;
  /* If we have two baryons in the beginning, we must have fermion resonance */
  if (type_particle1.pdgcode().baryon_number() != 0
      && type_particle2.pdgcode().baryon_number() != 0
      && !type_particle1.pdgcode().is_antiparticle_of(type_particle2.pdgcode())
      && type_resonance.pdgcode().baryon_number() == 0)
    return 0.0;

  /* Loop over particle types to find allowed combinations. */
  for (const ParticleType &second_type : ParticleType::list_all()) {
    /* We are interested only stable particles here */
    if (!second_type.is_stable()) {
      continue;
    }

    /* Check for charge conservation. */
    if (type_resonance.charge() + second_type.charge()
        != type_particle1.charge() + type_particle2.charge()) {
      continue;
    }

    /* Check for baryon number conservation. */
    int initial_baryon_number = type_particle1.pdgcode().baryon_number()
                              + type_particle2.pdgcode().baryon_number();
    int final_baryon_number = type_resonance.pdgcode().baryon_number()
                            + second_type.pdgcode().baryon_number();
    if (final_baryon_number != initial_baryon_number) {
      continue;
    }

    int I_z = type_resonance.isospin3() + second_type.isospin3();

    /* Compute total isospin range with given initial and final particles. */
    int I_max = std::min(type_resonance.isospin() + second_type.isospin(),
                         type_particle1.isospin() + type_particle2.isospin());
    int I_min = std::max(abs(type_resonance.isospin()
                             - second_type.isospin()),
                         abs(type_particle1.isospin()
                             - type_particle2.isospin()));
    I_min = std::max(I_min, abs(I_z));

    /* Loop over total isospin in allowed range.
     * Use decrement of 2, since isospin is multiplied by 2. */
    double isospin_factor = 0.;
    for (int I_tot = I_max; I_tot >= I_min; I_tot -= 2) {
      isospin_factor
        = isospin_clebsch_gordan(type_particle1, type_particle2, I_tot, I_z)
        * isospin_clebsch_gordan(type_resonance, second_type, I_tot, I_z);
    }

    /* If Clebsch-Gordan coefficient is zero, don't bother with the rest */
    if (std::abs(isospin_factor) < really_small) {
      continue;
    }

    /* Integration limits. */
    double lower_limit = type_resonance.minimum_mass();
    double upper_limit = std::sqrt(mandelstam_s) - second_type.mass();
    /* Check the available energy (requiring it to be a little above the
     * threshold, because the integration will not work if it's too close). */
    if (upper_limit - lower_limit < 1E-3) {
      continue;
    }

    /* Calculate resonance production cross section
     * using the Breit-Wigner distribution as probability amplitude.
     * Integrate over the allowed resonance mass range. */
    IntegrandParameters params = {&type_resonance, second_type.mass(),
                                  mandelstam_s};
    printd("Process: %s %s -> %s %s\n", type_particle1.name().c_str(),
     type_particle2.name().c_str(), second_type.name().c_str(),
     type_resonance.name().c_str());
    printd("Limits: %g %g \n", lower_limit, upper_limit);
    double resonance_integral, integral_error;
    quadrature_1d(&spectral_function_integrand, &params,
                  lower_limit, upper_limit,
                  &resonance_integral, &integral_error);
    printd("Integral value: %g Error: %g \n", resonance_integral,
      integral_error);

    /* Cross section for 2->2 process with resonance in final state.
     * Based on Eq. (51) in PhD thesis of J. Weil (Giessen U., 2013)
     * http://geb.uni-giessen.de/geb/volltexte/2013/10253/
     */
    float xsection
      = isospin_factor * isospin_factor
        / mandelstam_s
        / std::sqrt(cm_momentum_squared)
        * resonance_integral
        * nn_to_resonance_matrix_element(mandelstam_s, type_resonance,
                                         second_type);

    if (xsection > really_small) {
      process_list->push_back(ProcessBranch(type_resonance.pdgcode(),
                                            second_type.pdgcode(), xsection));
      number_of_processes++;
    }
  }
  return number_of_processes;
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

/* Scattering matrix amplitude squared.
 * Introduces energy-dependent modification
 * on the constant matrix element from
 * M. Effenberger diploma thesis (Giessen 1996)
 */
float nn_to_resonance_matrix_element(const double mandelstam_s,
  const ParticleType &type_final_a, const ParticleType &type_final_b) {
  PdgCode delta = PdgCode("2224");
  if (type_final_a.pdgcode().iso_multiplet()
      != type_final_b.pdgcode().iso_multiplet()) {
    /* N + N -> N + Delta: fit to Dmitriev OBE model,
       Nucl. Phys. A 459, 503 (1986) */
    if (type_final_a.pdgcode().iso_multiplet() == delta.iso_multiplet()
        || type_final_b.pdgcode().iso_multiplet() == delta.iso_multiplet()) {
      return 459. / std::pow(std::sqrt(mandelstam_s) - 1.104, 1.951);
    } else {
      return 0.0;
    }
  } else {
    return 0.0;
  }
}

}  // namespace Smash
