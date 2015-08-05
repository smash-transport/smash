/*
 *    Copyright (c) 2012-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#include <cmath>
#include <iostream>

#include "include/fourvector.h"
#include "include/iomanipulators.h"
#include "include/numerics.h"

namespace Smash {

/**
 * Algorithmic
 * -----------
 *
 * (Notation: \f$\vec{a}\f$ is a Three-Vector, \f$a^\mu\f$ is a Four-Vector.)
 *
 * The gamma factor is \f$\gamma = 1/\sqrt{1-\vec{v}^2}\f$.
 *
 * The time-like component of a Lorentz-boosted FourVector \f$x^\mu =
 * (x_0, x_1, x_2, x_3) = (x_0, \vec{r})\f$ with velocity \f$\vec v\f$
 * is
 *
 * \f{eqnarray*}{x^\prime_0 = \gamma \cdot (x_0 - \vec{r}\cdot\vec{v}),\f}
 *
 * and the space-like components i = 1, 2, 3 are:
 * \f{eqnarray*}{
 * x^\prime_i &=& x_i + v_i \cdot (\frac{\gamma - 1}{\vec{v}^2} \cdot \vec{r}\cdot\vec{v} - \gamma \cdot x_0)\\
 *            &=& x_i + v_i \cdot (\frac{\gamma^2}{\gamma + 1} \cdot \vec{r}\cdot\vec{v} - \gamma \cdot x_0)\\
 *            &=& x_i - v_i \cdot \frac{\gamma}{\gamma + 1} \cdot (\gamma(x_0 - \vec{r}\cdot\vec{v}) + x_0 )\\
 *            &=& x_i - v_i \cdot \frac{\gamma}{\gamma + 1} \cdot (x^\prime_0 + x_0)
 * \f}
 *
 * Note: This function is equivalent to -velocity Boost from ROOT
 */
FourVector FourVector::LorentzBoost(const ThreeVector &v) const {
  const double velocity_squared = v.sqr();

  const double gamma = velocity_squared < 1. ?
                       1. / std::sqrt(1. - velocity_squared) : 0;

  // this is used four times in the Vector:
  const double xprime_0 = gamma * (this->x0() - this->threevec()*v);
  // this is the part of the space-like components that is always the same:
  const double constantpart = gamma / (gamma + 1) * (xprime_0 + this->x0());
  return FourVector (xprime_0, this->threevec() - v * constantpart);
}

/* Check if all four vector components are almost equal
 * (accuracy \f$10^{-12}\f$). */
bool FourVector::operator==(const FourVector &a) const {
  return almost_equal(x_[0], a.x_[0])
      && almost_equal(x_[1], a.x_[1])
      && almost_equal(x_[2], a.x_[2])
      && almost_equal(x_[3], a.x_[3]);
}

std::ostream& operator<<(std::ostream& out, const FourVector& vec) {
  out.put('(');
  out.fill(' ');
  for (auto x : vec) {
    out << field<8> << x;
  }
  return out << ')';
}

}  // namespace Smash
