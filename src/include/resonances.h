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
 * **NB: This function assumes stable decay products!**
 */
float calculate_minimum_mass(Particles *particles, int pdgcode);

/**
 * Find all resonances and their production cross sections.
 *
 * Given two particles, create a list of possible resonance
 * production processes and their cross sections.
 * Process can be either 2-to-1 or 2-to-2.
 *
 * \return A list of processes with resonance in the final state.
 *  Each element in the list contains the type(s)
 *  of the final state particle(s)
 *  and the cross section for that particular process.
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
 *
 * \return The number of possible processes. Also adds elements
 * to the process_list.
 */
size_t two_to_two_formation(Particles *particles,
  const ParticleType &type_particle1,
  const ParticleType &type_particle2, const ParticleType &type_resonance,
  double mandelstam_s, double cm_momentum_squared,
  std::vector<ProcessBranch> *process_list);

/**
 * Function for 1-dimensional GSL integration.
 *
 * \param integrand_function Function of 1 variable to be integrated over.
 * \param parameters Container for possible parameters needed by the integrand.
 * \param lower_limit Lower limit of the integral.
 * \param upper_limit Upper limit of the integral.
 * \param integral_value Result of integration.
 * \param integral_error Uncertainty of the result.
 */
void quadrature_1d(double (*integrand_function)(double, void*),
                   std::vector<double> *parameters,
                   double lower_limit, double upper_limit,
                   double *integral_value, double *integral_error);

/// Spectral function of the resonance.
double spectral_function(double resonance_mass, double resonance_pole,
                         double resonance_width);

/**
 * Spectral function integrand for GSL integration.
 *
 * The integrand is \f$2m A(m) p_{cm}^f\f$, where \f$m\f$ is the
 * resonance mass, \f$A(m)\f$ is the spectral function
 *  and \f$p_{cm}^f\f$ is the center-of-mass momentum of the final state.
 *
 * \param resonance_mass Actual mass of the resonance.
 * \param parameters Container for the parameters needed
 * by the spectral function: Width of the resonance,
 * pole mass of the resonance, mass of the stable particle in the final state
 * and mandelstam-s of the process.
 */
double spectral_function_integrand(double resonance_mass, void * parameters);

/**
 * Resonance mass sampling for 2-particle final state
 * with *one resonance* and one *stable* particle.
 */
double sample_resonance_mass(Particles *particles, int pdg_resonance,
  int pdg_stable, double cms_energy);

/**
 * Resonance formation process.
 *
 * Creates one or two new particles, of which
 * one is a resonance.
 *
 * \param particles Particles in the simulation.
 * \param particle_id ID of the first initial state particle.
 * \param other_id ID of the second initial state particle.
 * \param produced_particles Final state particle type(s).
 * \return ID of the (first) new particle.
 */
int resonance_formation(Particles *particles, int particle_id, int other_id,
  std::vector<int> produced_particles);

}  // namespace Smash

#endif  // SRC_INCLUDE_RESONANCES_H_
