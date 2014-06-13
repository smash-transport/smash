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
 * Get the mass-dependent total width of particle with type t and mass m.
 * 
 * \param t Type of the decaying particle.
 * \param m Invariant mass of the decaying particle.
 */
float width_total(const ParticleType *t, const float m);

}  // namespace Smash

#endif  // SRC_INCLUDE_WIDTH_H_
