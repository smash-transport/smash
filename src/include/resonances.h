/*
 *    Copyright (c) 2013-2014
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

#include "forwarddeclarations.h"

#include <cstddef>
#include <cstdint>

namespace Smash {


/**
 * Calculate isospin Clebsch-Gordan coefficient
 *
 * \f$(-1)^{j_1 - j_2 + m_3} \sqrt(2 j_3 + 1) \cdot [Wigner 3J symbol] \f$
 * Note that the calculation assumes that isospin values
 * have been multiplied by two
 */
double clebsch_gordan_coefficient(const int isospin_a,
  const int isospin_b, const int isospin_resonance,
  const int isospin_z_a, const int isospin_z_b,
  const int isospin_z_resonance);


/**
 * Find all resonances that can be produced in a collision of the two
 * input particles, and the production cross sections of these resonances.
 *
 * Given the data and type information of two colliding particles,
 * create a list of possible resonance production processes
 * and their cross sections.
 * Process can be either 2-to-1 (just a resonance in the final state)
 * or 2-to-2 (resonance and a stable particle in the final state).
 *
 * \param[in] particle1 Data of the first colliding particle.
 * \param[in] particle2 Data of the second colliding particle.
 *
 * \return A list of processes with resonance in the final state.
 * Each element in the list contains the type(s)
 * of the final state particle(s)
 * and the cross section for that particular process.
 */
ProcessBranchList resonance_cross_section(const ParticleData &particle1,
                                          const ParticleData &particle2);

/**
 * Given the types of the two initial particles and a resonance,
 * return the 2-to-1 resonance production cross section.
 *
 * Checks are processed in the following order:
 * 1. Charge conservation
 * 2. Baryon number conservation
 * 3. Clebsch-Gordan
 * 4. Enough energy for all decay channels to be available for the resonance
 * 5. Detailed balance (reverse process exists)
 *
 * \param[in] type_particle1 Type information for the first initial particle.
 * \param[in] type_particle2 Type information for the second initial particle.
 * \param[in] type_resonance Type information for the resonance to be produced.
 * \param[in] mandelstam_s Mandelstam-s of the collision
 * of the two initial particles.
 * \param[in] cm_momentum_squared Square of the center-of-mass momentum of the
 * two initial particles.
 *
 * \return The cross section for the process
 * [initial particle 1] + [initial particle 2] -> resonance.
 */
double two_to_one_formation(const ParticleType &type_particle1,
                            const ParticleType &type_particle2,
                            const ParticleType &type_resonance,
                            double mandelstam_s, double cm_momentum_squared);

/**
 * Given the types of the two initial particles and a resonance,
 * return the final states containing the resonance and a stable particle
 * which are possible in the collision of the initial particles.
 *
 * Checks are processed in the following order:
 * 1. Charge conservation
 * 2. Baryon number conservation
 * 3. Clebsch-Gordan
 * 4. Enough energy for all decay channels to be available for the resonance
 *
 * \param[in] type_particle1 Type information for the first initial particle.
 * \param[in] type_particle2 Type information for the second initial particle.
 * \param[in] type_resonance Type information for the resonance to be produced.
 * \param[in] mandelstam_s Mandelstam-s of the collision
 * of the two initial particles.
 * \param[in] cm_momentum_squared Square of the center-of-mass momentum of the
 * two initial particles.
 * \param[in,out] process_list List of resonance production processes possible
 * in the collision of the two initial particles.
 * Each element in the list contains the type(s)
 * of the final state particle(s)
 * and the cross section for that particular process.
 *
 * \return The number of possible processes. Also adds elements
 * to the process_list.
 */
std::size_t two_to_two_formation(const ParticleType &type_particle1,
                                 const ParticleType &type_particle2,
                                 const ParticleType &type_resonance,
                                 double mandelstam_s,
                                 double cm_momentum_squared,
                                 ProcessBranchList *process_list);

/**
 * Spectral function
 * \f$A(m)=\frac{1}{\pi}\frac{m\Gamma(m)}{(m^2-m_0^2)^2+(m\Gamma(m))^2}\f$
 * of the resonance.
 */
double spectral_function(double resonance_mass, double resonance_pole,
                         double resonance_width);

/**
 * Spectral function integrand for GSL integration.
 *
 * The integrand is \f$2m A(m) p_{cm}^f\f$, where \f$m\f$ is the
 * resonance mass, \f$A(m)\f$ is the spectral function
 *  and \f$p_{cm}^f\f$ is the center-of-mass momentum of the final state.
 *
 * \param[in] resonance_mass Actual mass of the resonance.
 * \param[in] parameters Container for the parameters needed
 * by the spectral function: Width of the resonance,
 * pole mass of the resonance, mass of the stable particle in the final state
 * and mandelstam-s of the process.
 */
double spectral_function_integrand(double resonance_mass, void * parameters);

/**
 * Resonance mass sampling for 2-particle final state
 * with *one resonance* and one *stable* particle.
 *
 * \param[in] pdg_resonance PDG code of the resonance particle.
 * \param[in] pdg_stable PDG code of the stable particle.
 * \param[in] cms_energy center-of-mass energy
 * of the 2-particle final state.
 *
 * \return The mass of the resonance particle.
 */
double sample_resonance_mass(const ParticleType &type_resonance,
                             const ParticleType &type_stable,
                             const float cms_energy);

/**
 * Scattering matrix amplitude squared for \f$NN \rightarrow \f$
 *  \[resonant state\] processes.
 *
 * \param[in] mandelstam_s Mandelstam-s, i.e. collision CMS energy squared.
 * \param[in] type_final_a Type information for the first final state particle.
 * \param[in] type_final_b Type information for the second final state particle.
 *
 * \return Matrix amplitude squared \f$|\mathcal{M}(\sqrt{s})|^2\f$ for
 * process \f$N+N \rightarrow\f$ type_final_a + type_final_b
 */
float nn_to_resonance_matrix_element(const double mandelstam_s,
  const ParticleType &type_final_a, const ParticleType &type_final_b);

}  // namespace Smash

#endif  // SRC_INCLUDE_RESONANCES_H_
