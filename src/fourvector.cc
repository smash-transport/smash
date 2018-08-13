/*
 *    Copyright (c) 2012-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#include <cmath>
#include <iostream>

#include "smash/fourvector.h"
#include "smash/iomanipulators.h"
#include "smash/numerics.h"

namespace smash {

FourVector FourVector::LorentzBoost(const ThreeVector& v) const {
  const double velocity_squared = v.sqr();

  const double gamma =
      velocity_squared < 1. ? 1. / std::sqrt(1. - velocity_squared) : 0;

  // this is used four times in the Vector:
  const double xprime_0 = gamma * (this->x0() - this->threevec() * v);
  // this is the part of the space-like components that is always the same:
  const double constantpart = gamma / (gamma + 1) * (xprime_0 + this->x0());
  return FourVector(xprime_0, this->threevec() - v * constantpart);
}

bool FourVector::operator==(const FourVector& a) const {
  return almost_equal_physics(x_[0], a.x_[0]) &&
         almost_equal_physics(x_[1], a.x_[1]) &&
         almost_equal_physics(x_[2], a.x_[2]) &&
         almost_equal_physics(x_[3], a.x_[3]);
}

std::ostream& operator<<(std::ostream& out, const FourVector& vec) {
  out.put('(');
  out.fill(' ');
  for (auto x : vec) {
    out << field<8> << x;
  }
  return out << ')';
}

}  // namespace smash
