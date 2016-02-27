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

#include "particletype.h"

namespace Smash {

/**
 * Spectral function integrand for GSL integration, with one resonance in the
 * final state (the second particle is stable).
 *
 * The integrand is \f$ A(m) p_{cm}^f \f$, where \f$ m \f$ is the
 * resonance mass, \f$ A(m) \f$ is the spectral function
 *  and \f$ p_{cm}^f \f$ is the center-of-mass momentum of the final state.
 *
 * \param[in] resonance_mass Actual mass of the resonance.
 * \param[in] sqrts Center-of-mass energy, i.e. sqrt of Mandelstam s.
 * \param[in] stable_mass mass of the stable particle in the final state
 * \param[in] type type of the resonance
 */
float spec_func_integrand_1res(float resonance_mass, float sqrts,
                               float stable_mass, const ParticleType &type);


/**
 * Spectral function integrand for GSL integration, with two resonances in the
 * final state.
 *
 * The integrand is \f$ A_1(m_1) A_2(m_2) p_{cm}^f \f$, where \f$ m_1 \f$ and
 * \f$ m_2 \f$ are the resonance masses, \f$ A_1 \f$ and \f$ A_2 \f$ are the
 * spectral functions and \f$ p_{cm}^f \f$ is the center-of-mass momentum of
 * the final state.
 *
 * \param[in] sqrts Center-of-mass energy, i.e. sqrt of Mandelstam s.
 * \param[in] res_mass_1 Actual mass of the first resonance.
 * \param[in] res_mass_2 Actual mass of the second resonance.
 * \param[in] t1 Type of the first resonance.
 * \param[in] t2 Type of the second resonance.
 */
float spec_func_integrand_2res(float sqrts, float res_mass_1, float res_mass_2,
                               const ParticleType &t1, const ParticleType &t2);


/**
 * Resonance mass sampling for 2-particle final state
 * with *one resonance* and one *stable* particle.
 *
 * \param[in] type_res Type of the resonance particle.
 * \param[in] mass_stable Mass of the stable particle.
 * \param[in] cms_energy center-of-mass energy of the 2-particle final state.
 * \param[in] L relative angular momentum of the final-state particles
 *
 * \return The mass of the resonance particle.
 */
float sample_resonance_mass(ParticleType &type_res, const float mass_stable,
                            const float cms_energy, int L = 0);

/**
 * Resonance mass sampling for 2-particle final state with two resonances.
 *
 * \param[in] type_res_1 Type of the first resonance.
 * \param[in] type_res_2 Type of the second resonance.
 * \param[in] cms_energy center-of-mass energy of the 2-particle final state.
 * \param[in] L relative angular momentum of the final-state particles
 *
 * \return The masses of the resonance particles.
 */
std::pair<float, float> sample_resonance_masses(const ParticleType &type_res_1,
                                                const ParticleType &type_res_2,
                                                const float cms_energy,
                                                int L = 0);


}  // namespace Smash

#endif  // SRC_INCLUDE_RESONANCES_H_
