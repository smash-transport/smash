/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/particledata.h"

#include "include/width.h"
#include <iostream>

namespace Smash {

std::ostream &operator<<(std::ostream &s, const ParticleData &p) {
  return s << "ParticleData{ id: " << p.id() << ", '" << p.type().name()
           << "', pos: " << p.position() << ", mom: " << p.momentum() << " }";
}

}  // namespace Smash
