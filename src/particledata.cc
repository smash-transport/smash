/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/particledata.h"

#include "include/particles.h"
#include "include/particletype.h"

namespace Smash {

const ParticleType &ParticleData::type(const Particles &particles) const {
  return particles.particle_type(pdgcode_);
}

}  // namespace Smash
