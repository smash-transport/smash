/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/threevector.h"

#include "include/iomanipulators.h"

namespace Smash {

std::ostream &operator<<(std::ostream &out, const ThreeVector &v) {
  out.put('(');
  out.fill(' ');
  for (auto x : v) {
    out << field<8> << x;
  }
  return out << ')';
}

}  // namespace Smash
