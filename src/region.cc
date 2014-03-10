/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/region.h"

namespace Smash {

Region::Region(const ParticleList &p) : particles_(p) {}

void Region::find_all_actions() {
  for (const auto &finder : all_action_finders_) {
    all_actions_ += finder.find_possible_actions(particles_);
  }
}

}  // namespace Smash
