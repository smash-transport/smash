/*
 *
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/particles.h"

namespace Smash {

/* boost_CM - boost to center of momentum */
ThreeVector boost_CM(ParticleData *particle1, ParticleData *particle2) {
  // determine CMS velocity
  FourVector mom = particle1->momentum() + particle2->momentum();
  ThreeVector velocity = mom.threevec() / mom.x0();

  // Boost the momenta into CMS frame
  particle1->boost(velocity);
  particle2->boost(velocity);

  return velocity;
}

/* boost_back_CM - boost back from center of momentum */
void boost_back_CM(ParticleData *particle1, ParticleData *particle2,
                   const ThreeVector &velocity_orig) {

  ThreeVector velocity = - velocity_orig;

  particle1->boost(velocity);
  particle2->boost(velocity);
}


Particles::Particles() {}


void Particles::reset() {
  id_max_ = -1;
  data_.clear();
}

}  // namespace Smash
