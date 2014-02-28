/*
 *
 *    Copyright (c) 2013
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

namespace Smash {

/* calculate_minimum_mass
 * - calculate the minimum rest energy the resonance must have
 * to be able to decay through any of its decay channels
 * NB: This function assumes stable decay products!
 */
float calculate_minimum_mass(Particles *particles, int pdgcode) {
  /* If the particle happens to be stable, just return the mass */
  if (particles->particle_type(pdgcode).width() < 0.0)
    return particles->particle_type(pdgcode).mass();
  /* Otherwise, let's find the highest mass value needed in any decay mode */
  float minimum_mass = 0.0;
  const std::vector<ProcessBranch> decaymodes
    = particles->decay_modes(pdgcode).decay_mode_list();
  for (std::vector<ProcessBranch>::const_iterator mode = decaymodes.begin();
       mode != decaymodes.end(); ++mode) {
    size_t decay_particles = mode->particle_list().size();
    float total_mass = 0.0;
    for (size_t i = 0; i < decay_particles; i++) {
      /* Stable decay products assumed; for resonances the mass can be lower! */
      total_mass
        += particles->particle_type(mode->particle_list().at(i)).mass();
    }
    if (total_mass > minimum_mass)
      minimum_mass = total_mass;
  }
  return minimum_mass;
}

/* resonance_cross_section - energy-dependent cross section
 * for producing a resonance
 */
std::vector<ProcessBranch> resonance_cross_section(
    const ParticleData &particle1, const ParticleData &particle2,
    const ParticleType &type_particle1, const ParticleType &type_particle2,
    Particles *particles) {
  std::vector<ProcessBranch> resonance_process_list;

  /* first item refers to total resonance production cross section */
  ProcessBranch resonance_process;
  resonance_process.add_particle(0);
  resonance_process.set_weight(0.0);
  resonance_process_list.push_back(resonance_process);

  /* Isospin symmetry factor, by default 1 */
  int symmetryfactor = 1;
  /* Do the particles have the same isospin value? */
  if (type_particle1.isospin() == type_particle2.isospin()) {
    /* Do they have the same spin? */
    if (type_particle1.spin() == type_particle2.spin()) {
      /* Are their PDG codes of same length? */
      int abs_pdg1 = abs(type_particle1.pdgcode()), digits1 = 0;
      while (abs_pdg1) {
        abs_pdg1 /= 10;
        digits1++;
      }
      int abs_pdg2 = abs(type_particle2.pdgcode()), digits2 = 0;
      while (abs_pdg2) {
        abs_pdg2 /= 10;
        digits2++;
      }
      if (digits1 == digits2) {
        /* If baryons, do they have the same baryon number? */
        if (type_particle1.spin() % 2 == 0 ||
            std::signbit(type_particle1.pdgcode())
            == std::signbit(type_particle2.pdgcode())) {
          /* Ok, particles are in the same isospin multiplet,
             apply symmetry factor */
          symmetryfactor = 2;
        }
      }
    }
  }

  /* Mandelstam s = (p_a + p_b)^2 = square of CMS energy */
  const double mandelstam_s =
       (particle1.momentum() + particle2.momentum()).Dot(
         particle1.momentum() + particle2.momentum());

  /* CM momentum */
  const double cm_momentum_squared
    = (particle1.momentum().Dot(particle2.momentum())
       * particle1.momentum().Dot(particle2.momentum())
       - type_particle1.mass() * type_particle1.mass()
       * type_particle2.mass() * type_particle2.mass()) / mandelstam_s;

  /* Find all the possible resonances */
  for (auto i = particles->types_cbegin(); i != particles->types_cend(); ++i) {
    ParticleType type_resonance = i->second;

    /* Not a resonance, go to next type of particle */
    if (type_resonance.width() < 0.0)
      continue;

    /* Same resonance as in the beginning, ignore */
    if ((type_particle1.width() > 0.0
         && type_resonance.pdgcode() == type_particle1.pdgcode())
        || (type_particle2.width() > 0.0
            && type_resonance.pdgcode() == type_particle2.pdgcode()))
      continue;

    /* No decay channels found, ignore */
    if (particles->decay_modes(type_resonance.pdgcode()).empty())
      continue;

    float resonance_xsection
      = symmetryfactor * two_to_one_formation(particles, type_particle1,
        type_particle2, type_resonance, mandelstam_s, cm_momentum_squared);

    /* If cross section is non-negligible, add resonance to the list */
    if (resonance_xsection > really_small) {
      resonance_process.clear();
      resonance_process.add_particle(type_resonance.pdgcode());
      resonance_process.set_weight(resonance_xsection);
      resonance_process_list.push_back(resonance_process);
      resonance_process_list.at(0).change_weight(resonance_xsection);

      printd("Found resonance %i (%s) with mass %f and width %f.\n",
             type_resonance.pdgcode(), type_resonance.name().c_str(),
             type_resonance.mass(), type_resonance.width());
      printd("2->1 with original particles: %s %s Charges: %i %i \n",
             type_particle1.name().c_str(), type_particle2.name().c_str(),
             type_particle1.charge(), type_particle2.charge());
    }
    /* Same procedure for possible 2->2 resonance formation processes */
    /* XXX: For now, we allow this only for baryon-baryon interactions */
    if (type_particle1.spin() % 2 != 0 && type_particle2.spin() % 2 != 0) {
      size_t two_to_two_processes
         = two_to_two_formation(particles, type_particle1,
           type_particle2, type_resonance, mandelstam_s, cm_momentum_squared,
           symmetryfactor, &resonance_process_list);
      if (two_to_two_processes > 0) {
        printd("Found %zu 2->2 processes for resonance %i (%s).\n",
               two_to_two_processes,
               type_resonance.pdgcode(), type_resonance.name().c_str());
        printd("2->2 with original particles: %s %s Charges: %i %i \n",
               type_particle1.name().c_str(), type_particle2.name().c_str(),
               type_particle1.charge(), type_particle2.charge());
      }
    }
  }
  return resonance_process_list;
}

/* two_to_one_formation -- only the resonance in the final state */
double two_to_one_formation(Particles *particles,
  const ParticleType &type_particle1,
  const ParticleType &type_particle2, const ParticleType &type_resonance,
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
    if (type_particle1.spin() % 2 != 0
        && (std::signbit(type_particle1.pdgcode())
            != std::signbit(type_resonance.pdgcode()))) {
      return 0.0;
    } else if (type_particle2.spin() % 2 != 0
        && (std::signbit(type_particle2.pdgcode())
        != std::signbit(type_resonance.pdgcode()))) {
      return 0.0;
    }
  }

  /* Isospin z-component based on Gell-Mann–Nishijima formula
   * 2 * Iz = 2 * charge - (baryon number + strangeness + charm)
   * XXX: Strangeness and charm ignored for now!
   */
  const int isospin_z1 = type_particle1.spin() % 2 == 0
    ? type_particle1.charge() * 2
    : type_particle1.charge() * 2 - type_particle1.pdgcode()
                                    / abs(type_particle1.pdgcode());
  const int isospin_z2 = type_particle2.spin() % 2 == 0
    ? type_particle2.charge() * 2
    : type_particle2.charge() * 2 - type_particle2.pdgcode()
                                    / abs(type_particle2.pdgcode());
  int isospin_z_resonance = (type_resonance.spin()) % 2 == 0
    ? type_resonance.charge() * 2
    : type_resonance.charge() * 2 - type_resonance.pdgcode()
                                    / abs(type_resonance.pdgcode());

  /* Calculate isospin Clebsch-Gordan coefficient
   * (-1)^(j1 - j2 + m3) * sqrt(2 * j3 + 1) * [Wigner 3J symbol]
   * Note that the calculation assumes that isospin values
   * have been multiplied by two
   */
  double wigner_3j =  gsl_sf_coupling_3j(type_particle1.isospin(),
     type_particle2.isospin(), type_resonance.isospin(),
     isospin_z1, isospin_z2, -isospin_z_resonance);
  double clebsch_gordan_isospin = 0.0;
  if (fabs(wigner_3j) > really_small)
    clebsch_gordan_isospin = pow(-1, type_particle1.isospin() / 2.0
    - type_particle2.isospin() / 2.0 + isospin_z_resonance / 2.0)
    * sqrt(type_resonance.isospin() + 1) * wigner_3j;

  printd("CG: %g I1: %i I2: %i IR: %i iz1: %i iz2: %i izR: %i \n",
       clebsch_gordan_isospin,
       type_particle1.isospin(), type_particle2.isospin(),
       type_resonance.isospin(),
       isospin_z1, isospin_z2, isospin_z_resonance);

  /* If Clebsch-Gordan coefficient is zero, don't bother with the rest */
  if (fabs(clebsch_gordan_isospin) < really_small)
    return 0.0;

  /* Check the decay modes of this resonance */
  const std::vector<ProcessBranch> decaymodes
    = particles->decay_modes(type_resonance.pdgcode()).decay_mode_list();
  bool not_enough_energy = false;
  /* Detailed balance required: Formation only possible if
   * the resonance can decay back to these particles
   */
  bool not_balanced = true;
  for (std::vector<ProcessBranch>::const_iterator mode
       = decaymodes.begin(); mode != decaymodes.end(); ++mode) {
    size_t decay_particles = mode->particle_list().size();
    if ( decay_particles > 3 ) {
      printf("Warning: Not a 1->2 or 1->3 process!\n");
      printf("Number of decay particles: %zu \n", decay_particles);
    } else {
      /* There must be enough energy to produce all decay products */
      float mass_a, mass_b, mass_c = 0.0;
      mass_a = calculate_minimum_mass(particles, mode->particle_list().at(0));
      mass_b = calculate_minimum_mass(particles, mode->particle_list().at(1));
      if (decay_particles == 3) {
        mass_c = calculate_minimum_mass(particles, mode->particle_list().at(2));
      }
      if (sqrt(mandelstam_s) < mass_a + mass_b + mass_c)
        not_enough_energy = true;
      /* Initial state is also a possible final state;
       * weigh the cross section with the ratio of this branch
       * XXX: For now, assuming only 2-particle initial states
       */
      if (decay_particles == 2
          && ((mode->particle_list().at(0) == type_particle1.pdgcode()
               && mode->particle_list().at(1) == type_particle2.pdgcode())
              || (mode->particle_list().at(0) == type_particle2.pdgcode()
                  && mode->particle_list().at(1) == type_particle1.pdgcode()))
          && (mode->weight() > 0.0))
        not_balanced = false;
    }
  }
  if (not_enough_energy)
    return 0.0;

  if (not_balanced)
    return 0.0;

  /* Calculate spin factor */
  const double spinfactor = (type_resonance.spin() + 1)
    / ((type_particle1.spin() + 1) * (type_particle2.spin() + 1));
  float resonance_width = type_resonance.width();
  float resonance_mass = type_resonance.mass();
  /* Calculate resonance production cross section
   * using the Breit-Wigner distribution as probability amplitude
   * See Eq. (176) in Buss et al., Physics Reports 512, 1 (2012)
   */
  return clebsch_gordan_isospin * clebsch_gordan_isospin * spinfactor
         * 4.0 * M_PI / cm_momentum_squared
         * breit_wigner(mandelstam_s, resonance_mass, resonance_width)
         * hbarc * hbarc / fm2_mb;
}

/* two_to_two_formation -- resonance and another particle in final state */
size_t two_to_two_formation(Particles *particles,
  const ParticleType &type_particle1,
  const ParticleType &type_particle2, const ParticleType &type_resonance,
  double mandelstam_s, double cm_momentum_squared, double symmetryfactor,
  std::vector<ProcessBranch> *process_list) {
  size_t number_of_processes = 0;
  /* If we have two baryons in the beginning, we must have fermion resonance */
  if (type_particle1.spin() % 2 != 0 && type_particle2.spin() % 2 != 0
      && type_particle1.pdgcode() != -type_particle2.pdgcode()
      && type_resonance.spin() % 2 == 0)
    return 0.0;

  /* Isospin z-component based on Gell-Mann–Nishijima formula
   * 2 * Iz = 2 * charge - (baryon number + strangeness + charm)
   * XXX: Strangeness and charm ignored for now!
   */
  const int isospin_z1 = type_particle1.spin() % 2 == 0
    ? type_particle1.charge() * 2
    : type_particle1.charge() * 2 - type_particle1.pdgcode()
                                    / abs(type_particle1.pdgcode());
  const int isospin_z2 = type_particle2.spin() % 2 == 0
    ? type_particle2.charge() * 2
    : type_particle2.charge() * 2 - type_particle2.pdgcode()
                                    / abs(type_particle2.pdgcode());

  int isospin_z_resonance = (type_resonance.spin()) % 2 == 0
    ? type_resonance.charge() * 2
    : type_resonance.charge() * 2 - type_resonance.pdgcode()
                                    / abs(type_resonance.pdgcode());

  /* Compute initial total combined isospin range */
  int initial_total_maximum
    = type_particle1.isospin() + type_particle2.isospin();
  int initial_total_minimum
    = abs(type_particle1.isospin() - type_particle2.isospin());
  /* Loop over particle types to find allowed combinations */
  for (std::map<int, ParticleType>::const_iterator
       type_i = particles->types_cbegin(); type_i != particles->types_cend();
        ++type_i) {
    /* We are interested only stable particles here */
    if (type_i->second.width() > 0.0)
      continue;

    /* Check for charge conservation */
    if (type_resonance.charge() + type_i->second.charge()
        != type_particle1.charge() + type_particle2.charge())
      continue;

    /* Check for baryon number conservation */
    int initial_baryon_number = 0;
    if (type_particle1.spin() % 2 != 0) {
      initial_baryon_number += type_particle1.pdgcode()
                               / abs(type_particle1.pdgcode());
    }
    if (type_particle2.spin() % 2 != 0) {
      initial_baryon_number += type_particle2.pdgcode()
                               / abs(type_particle2.pdgcode());
    }
    int final_baryon_number = 0;
    if (type_resonance.spin() % 2 != 0) {
      final_baryon_number += type_resonance.pdgcode()
                               / abs(type_resonance.pdgcode());
    }
    if (type_i->second.spin() % 2 != 0) {
      final_baryon_number += type_i->second.pdgcode()
                               / abs(type_i->second.pdgcode());
    }
    if (final_baryon_number != initial_baryon_number)
      continue;

    /* Compute total isospin range with given initial and final particles */
    int isospin_maximum = std::min(type_resonance.isospin()
      + type_i->second.isospin(), initial_total_maximum);
    int isospin_minimum = std::max(abs(type_resonance.isospin()
      - type_i->second.isospin()), initial_total_minimum);

    int isospin_z_i = (type_i->second.spin()) % 2 == 0
    ? type_i->second.charge() * 2
    : type_i->second.charge() * 2 - type_i->second.pdgcode()
       / abs(type_i->second.pdgcode());
    int isospin_z_final = isospin_z_resonance + isospin_z_i;

    int isospin_final = isospin_maximum;
    double clebsch_gordan_isospin = 0.0;
    while (isospin_final >= isospin_minimum) {
      if (abs(isospin_z_final) > isospin_final)
        break;
      /* Calculate isospin Clebsch-Gordan coefficient combinations
       * (-1)^(j1 - j2 + m3) * sqrt(2 * j3 + 1) * [Wigner 3J symbol]
       * Note that the calculation assumes that isospin values
       * have been multiplied by two
       */
      double wigner_3j =  gsl_sf_coupling_3j(type_particle1.isospin(),
        type_particle2.isospin(), isospin_final,
        isospin_z1, isospin_z2, -isospin_z_final);
      if (fabs(wigner_3j) > really_small)
        clebsch_gordan_isospin += pow(-1, type_particle1.isospin() / 2.0
        - type_particle2.isospin() / 2.0 + isospin_z_final / 2.0)
        * sqrt(isospin_final + 1) * wigner_3j;

      printd("CG: %g I1: %i I2: %i IR: %i iz1: %i iz2: %i izR: %i \n",
         clebsch_gordan_isospin,
         type_particle1.isospin(), type_particle2.isospin(),
         isospin_final,
         isospin_z1, isospin_z2, isospin_z_final);
      /* isospin is multiplied by 2,
       *  so we must also decrease it by increments of 2
       */
      isospin_final = isospin_final - 2;
    }
    /* If Clebsch-Gordan coefficient is zero, don't bother with the rest */
    if (fabs(clebsch_gordan_isospin) < really_small)
      continue;

    /* Check the decay modes of this resonance */
    const std::vector<ProcessBranch> decaymodes
      = particles->decay_modes(type_resonance.pdgcode()).decay_mode_list();
    bool not_enough_energy = false;
    double minimum_mass = 0.0;
    for (std::vector<ProcessBranch >::const_iterator mode
         = decaymodes.begin(); mode != decaymodes.end(); ++mode) {
      size_t decay_particles = mode->particle_list().size();
      if ( decay_particles > 3 ) {
        printf("Warning: Not a 1->2 or 1->3 process!\n");
        printf("Number of decay particles: %zu \n", decay_particles);
      } else {
        /* There must be enough energy to produce all decay products */
        float mass_a, mass_b, mass_c = 0.0;
        mass_a = calculate_minimum_mass(particles, mode->particle_list().at(0));
        mass_b = calculate_minimum_mass(particles, mode->particle_list().at(1));
        if (decay_particles == 3) {
          mass_c = calculate_minimum_mass(particles,
                     mode->particle_list().at(2));
        }
        if (sqrt(mandelstam_s) < mass_a + mass_b + mass_c
                                 + type_i->second.mass()) {
          not_enough_energy = true;
        } else if (minimum_mass < mass_a + mass_b + mass_c) {
          minimum_mass = mass_a + mass_b + mass_c;
        }
      }
    }
    if (not_enough_energy)
      continue;

    /* Calculate resonance production cross section
     * using the Breit-Wigner distribution as probability amplitude
     * Integrate over the allowed resonance mass range
     */
    std::vector<double> integrand_parameters;
    integrand_parameters.push_back(type_resonance.width());
    integrand_parameters.push_back(type_resonance.mass());
    integrand_parameters.push_back(type_i->second.mass());
    integrand_parameters.push_back(mandelstam_s);
    double lower_limit = minimum_mass;
    double upper_limit = (sqrt(mandelstam_s) - type_i->second.mass());
    printd("Process: %s %s -> %s %s\n", type_particle1.name().c_str(),
     type_particle2.name().c_str(), type_i->second.name().c_str(),
     type_resonance.name().c_str());
    printd("Limits: %g %g \n", lower_limit, upper_limit);
    double resonance_integral, integral_error;
    quadrature_1d(&spectral_function_integrand, &integrand_parameters,
                  lower_limit, upper_limit,
                  &resonance_integral, &integral_error);
    printd("Integral value: %g Error: %g \n", resonance_integral,
      integral_error);

    /* matrix element squared over 16pi (in mb GeV^2)
     * (uniform angular distribution assumed)
     */
    double matrix_element = 180;

    /* Cross section for 2->2 process with resonance in final state.
     * Based on the general differential form in
     * Buss et al., Physics Reports 512, 1 (2012), Eq. (D.28)
     */
    float xsection = clebsch_gordan_isospin * clebsch_gordan_isospin
                      * symmetryfactor
                      * matrix_element
                      / mandelstam_s
                      / sqrt(cm_momentum_squared)
                      * resonance_integral;

    if (xsection > really_small) {
      ProcessBranch final_state;
      final_state.add_particle(type_resonance.pdgcode());
      final_state.add_particle(type_i->second.pdgcode());
      final_state.set_weight(xsection);
      process_list->push_back(final_state);
      process_list->at(0).change_weight(xsection);
      number_of_processes++;
    }
  }
  return number_of_processes;
}

/* Function for 1-dimensional GSL integration  */
void quadrature_1d(double (*integrand_function)(double, void*),
                     std::vector<double> *parameters,
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
  std::vector<double> *integrand_parameters
    = reinterpret_cast<std::vector<double> *>(parameters);
  double resonance_width = integrand_parameters->at(0);
  double resonance_pole_mass = integrand_parameters->at(1);
  double stable_mass = integrand_parameters->at(2);
  double mandelstam_s = integrand_parameters->at(3);

  /* center-of-mass momentum of final state particles */
  if (mandelstam_s - (stable_mass + resonance_mass)
      * (stable_mass + resonance_mass) > 0.0) {
    double cm_momentum_final
      = sqrt((mandelstam_s - (stable_mass + resonance_mass)
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
double sample_resonance_mass(Particles *particles, int pdg_resonance,
  int pdg_stable, double cms_energy) {
  /* First, find the minimum mass of this resonance */
  double minimum_mass
    = calculate_minimum_mass(particles, pdg_resonance);
  /* Define distribution parameters */
  float mass_stable
    = particles->particle_type(pdg_stable).mass();
  std::vector<double> parameters;
  parameters.push_back(particles->particle_type(pdg_resonance).width());
  parameters.push_back(particles->particle_type(pdg_resonance).mass());
  parameters.push_back(mass_stable);
  parameters.push_back(cms_energy * cms_energy);

  /* sample resonance mass from the distribution
   * used for calculating the cross section
   */
  double mass_resonance = 0.0, random_number = 1.0;
  double distribution_max
    = spectral_function_integrand(parameters.at(0), &parameters);
  double distribution_value = 0.0;
  while (random_number > distribution_value) {
    random_number = distribution_max * drand48();
    mass_resonance = (cms_energy - mass_stable - minimum_mass) * drand48()
                     + minimum_mass;
    distribution_value
      = spectral_function_integrand(mass_resonance, &parameters);
  }
  return mass_resonance;
}

/* Resonance formation kinematics */
int resonance_formation(Particles *particles, int particle_id, int other_id,
                        std::vector<int> produced_particles) {
  if (produced_particles.empty()) {
    printf("resonance_formation:\n");
    printf("Warning: No final state particles found!\n");
    printf("Resonance formation canceled. Returning -1.\n");
    return -1;
  }
  /* Template for a new particle */
  ParticleData resonance;

  const double cms_energy = particles->data(particle_id).momentum().x0()
    + particles->data(other_id).momentum().x0();

  int id_first_new = -1;
  if (produced_particles.size() == 1) {
    resonance.set_pdgcode(produced_particles.at(0));
    /* Center-of-momentum frame of initial particles
     * is the rest frame of the resonance
     *
     * We use fourvector to set 4-momentum, as setting it
     * with doubles requires that particle is on
     * mass shell, which is not generally true for resonances
     */
    FourVector resonance_momentum(cms_energy, 0.0, 0.0, 0.0);
    resonance.set_momentum(resonance_momentum);

    printd("Momentum of the new particle: %g %g %g %g \n",
      resonance.momentum().x0(),
      resonance.momentum().x1(),
      resonance.momentum().x2(),
      resonance.momentum().x3());

    /* Initialize position */
    resonance.set_position(1.0, 0.0, 0.0, 0.0);
    id_first_new = particles->add_data(resonance);
  } else if (produced_particles.size() == 2) {
    /* 2 particles in final state. Need another particle template */
    /* XXX: For now, it is assumed that the other particle is stable! */
    ParticleData stable_product;
    if (particles->particle_type(produced_particles.at(0)).width() > 0) {
      resonance.set_pdgcode(produced_particles.at(0));
      stable_product.set_pdgcode(produced_particles.at(1));
    } else {
      stable_product.set_pdgcode(produced_particles.at(0));
      resonance.set_pdgcode(produced_particles.at(1));
    }
    float mass_stable
      = particles->particle_type(stable_product.pdgcode()).mass();
    /* Sample resonance mass */
    double mass_resonance = sample_resonance_mass(particles,
      resonance.pdgcode(), stable_product.pdgcode(), cms_energy);

    /* Sample the particle momenta */
    sample_cms_momenta(&resonance, &stable_product, cms_energy, mass_resonance,
                       mass_stable);

    /* Initialize positions */
    resonance.set_position(1.0, 0.0, 0.0, 0.0);
    stable_product.set_position(1.0, 0.0, 0.0, 0.0);
    id_first_new = particles->add_data(resonance);
    particles->add_data(stable_product);
  } else {
    printf("resonance_formation:\n");
    printf("Warning: %zu particles in final state!\n",
           produced_particles.size());
    printf("Resonance formation canceled. Returning -1.\n");
    return -1;
  }
  /* Return the id of the first new particle */
  return id_first_new;
}

}  // namespace Smash
