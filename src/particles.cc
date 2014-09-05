/*
 *
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/particles.h"
#include <iostream>
#include <iomanip>

namespace Smash {

Particles::Particles() {}

void Particles::reset() {
  id_max_ = -1;
  data_.clear();
}

std::ostream &operator<<(std::ostream &out, const Particles &particles) {
  using namespace std;
  out << "Particles:\n";
  int n = 0;
  for (const auto &p : particles.data()) {
    while (p.id() != n) {
      out << "------  ";
      if ((++n & 15) == 0) {
        out << '\n';
      }
    }
    out << setw(5) << setprecision(3) << p.effective_mass() << p.type().name();
    if ((++n & 15) == 0) {
      out << '\n';
    }
  }
  return out;
}

}  // namespace Smash
