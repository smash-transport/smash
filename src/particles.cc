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

Particles::Particles() {}


void Particles::reset() {
  id_max_ = -1;
  data_.clear();
}

}  // namespace Smash
