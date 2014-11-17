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
 * \param mass_a Mass of the first daughter particle [GeV].
 * \param mass_b Mass of the second daughter particle [GeV].
 * \param L Angular momentum of the decay.
 * \param partialWidth_pole Partial width at the pole mass [GeV].
 */
float width_Manley_stable(const float mass, const float poleMass,
                          const float mass_a, const float mass_b,
                          const int L, const float partialWidth_pole);

/**
 * Get the mass-dependent width of a two-body decay into one stable and one
 * unstable particle according to Manley/Saleski, Phys. Rev. D 45 (1992) 4002.
 *
 * \param mass Actual mass of the decaying particle [GeV].
 * \param poleMass Pole mass of the decaying particle [GeV].
 * \param mass_stable Mass of the stable daughter particle [GeV].
 * \param type_unstable Type of the unstable daughter particle [GeV].
 * \param L Angular momentum of the decay.
 * \param partialWidth_pole Partial width at the pole mass [GeV].
 */
float width_Manley_semistable(const float mass, const float poleMass,
                              const float mass_stable,
                              const ParticleType *type_unstable,
                              const int L, const float partialWidth_pole);

}  // namespace Smash

#endif  // SRC_INCLUDE_WIDTH_H_
