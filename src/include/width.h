/*
 *    Copyright (c) 2012-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_WIDTH_H_
#define SRC_INCLUDE_WIDTH_H_

#include "particletype.h"

namespace Smash {

/**
 * Get the mass-dependent width of a two-body decay into stable particles
 * according to Manley/Saleski, Phys. Rev. D 45 (1992) 4002.
 * 
 * \param mass Actual mass of the decaying particle [GeV].
 * \param poleMass Pole mass of the decaying particle [GeV].
 * \param mass1 Mass of the first daughter particle [GeV].
 * \param mass2 Mass of the second daughter particle [GeV].
 * \param L Angular momentum of the decay.
 * \param partialWidth_pole Partial width at the pole mass [GeV].
 */
float width_Manley (const float mass, const float poleMass,
                    const float mass1, const float mass2,
                    const int L, const float partialWidth_pole);

float width_total (const ParticleType *t, const float m);

}  // namespace Smash

#endif  // SRC_INCLUDE_WIDTH_H_
