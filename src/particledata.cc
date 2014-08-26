/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/particledata.h"

#include "include/iomanipulators.h"
#include <iostream>

namespace Smash {

std::ostream &operator<<(std::ostream &out, const ParticleData &p) {
  using namespace std;
  out.fill(' ');
  return out << p.type().name() << right
             << "{id:" << field<6> << p.id()
             << ", pos [fm]:"  // XXX: is fm correct?
             << p.position() << ", mom [GeV]:" << p.momentum() << "}";
}

}  // namespace Smash
