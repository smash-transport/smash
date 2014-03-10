/*
 *    Copyright (c) 2013
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

/**
 * \file resonances.h
 * Functions related to resonance production.
 */

#ifndef SRC_INCLUDE_RESONANCES_H_
#define SRC_INCLUDE_RESONANCES_H_

#include <cstddef>
#include <vector>

namespace Smash {

/* necessary forward declarations */
class Particles;
class ParticleData;
class ParticleType;
class ProcessBranch;

/**
 * The minimum mass of the resonance.
 *
 * Calculate the minimum rest energy the resonance must have
 * to be able to decay through any of its decay channels.
 * NB: This function assumes stable decay products!
 */
float calculate_minimum_mass(Particles *particles, int pdgcode);

/**
 * Find all resonances and their production cross sections.
 *
 * Given two particles, create a list of possible resonance
 * production processes and their cross sections.
 * Process can be either 2-to-1 or 2-to-2.
 */
std::vector<ProcessBranch> resonance_cross_section(
  const ParticleData &particle1, const ParticleData &particle2,
  const ParticleType &type_particle1, const ParticleType &type_particle2,
  Particles *particles);

/**
 * Given two initial particles and a resonance,
 * compute the 2-to-1 resonance production cross section.
 */
double two_to_one_formation(Particles *particles,
  const ParticleType &type_particle1,
  const ParticleType &type_particle2, const ParticleType &type_resonance,
  double mandelstam_s, double cm_momentum_squared);

/**
 * Given two initial particles and a resonance,
 * compute cross sections for all possible 2-to-2 processes
 * with that resonance in the final state.
 */
size_t two_to_two_formation(Particles *particles,
  const ParticleType &type_particle1,
  const ParticleType &type_particle2, const ParticleType &type_resonance,
  double mandelstam_s, double cm_momentum_squared,
  std::vector<ProcessBranch> *process_list);

/// Function for 1-dimensional GSL integration.
void quadrature_1d(double (*integrand_function)(double, void*),
                   std::vector<double> *parameters,
                   double lower_limit, double upper_limit,
                   double *integral_value, double *integral_error);

/// Spectral function of the resonance.
double spectral_function(double resonance_mass, double resonance_pole,
                         double resonance_width);

/// Spectral function integrand for GSL integration.
double spectral_function_integrand(double mandelstam_s, void * parameters);

/**
 * Resonance mass sampling for 2-particle final state
 * with one resonance and one stable particle.
 */
double sample_resonance_mass(Particles *particles, int pdg_resonance,
  int pdg_stable, double cms_energy);

/// 2-to-1 resonance formation process.
int resonance_formation(Particles *particles, int particle_id, int other_id,
  std::vector<int> produced_particles);

}  // namespace Smash

#endif  // SRC_INCLUDE_RESONANCES_H_
