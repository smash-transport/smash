/*
 *
 *    Copyright (c) 2014-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/threevector.h"

#include "smash/iomanipulators.h"

namespace smash {

std::ostream &operator<<(std::ostream &out, const ThreeVector &v) {
  out.put('(');
  out.fill(' ');
  for (auto x : v) {
    out << field<8> << x;
  }
  return out << ')';
}

}  // namespace smash
