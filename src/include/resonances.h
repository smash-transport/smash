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
float sample_resonance_mass(const ParticleType &type_res,
                            const float mass_stable,
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
