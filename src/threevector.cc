/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/threevector.h"

#include <cmath>

#include "include/constants.h"
#include "include/iomanipulators.h"

namespace Smash {

double ThreeVector::abs() const {
  return std::sqrt((*this)*(*this));
}

double ThreeVector::get_phi() const {
  if (std::abs(x1()) < really_small && std::abs(x2()) < really_small) {
    return 0.;
  } else {
    return std::atan2(x2(), x1());
  }
}

double ThreeVector::get_theta() const {
  double r = abs();
  return (r > 0.) ? std::acos(x3()/r) : 0.;
}

void ThreeVector::rotate(double phi, double theta, double psi) {
  // Compute the cosine and sine for each angle.
  double cos_phi = std::cos(phi);
  double sin_phi = std::sin(phi);
  double cos_theta = std::cos(theta);
  double sin_theta = std::sin(theta);
  double cos_psi = std::cos(psi);
  double sin_psi = std::sin(psi);
  // Get original coordinates.
  std::array<double, 3> x_old = x_;
  // Compute new coordinates.
  x_[0] = (cos_phi * cos_psi - sin_phi * cos_theta * sin_psi) * x_old[0]
        + (sin_phi * cos_psi + cos_phi * cos_theta * sin_psi) * x_old[1]
        + sin_theta * sin_psi * x_old[2];
  x_[1] = (-cos_phi * sin_psi - sin_phi * cos_theta * cos_psi) * x_old[0]
        + (-sin_phi * sin_psi + cos_phi * cos_theta * cos_psi) * x_old[1]
        + sin_theta * cos_psi * x_old[2];
  x_[2] = sin_phi * sin_theta * x_old[0]
        - cos_phi * sin_theta * x_old[1]
        + cos_theta * x_old[2];
}

void ThreeVector::rotate_around_y(double theta) {
  double cost = std::cos(theta);
  double sint = std::sin(theta);
  // Get original coordinates.
  std::array<double, 3> x_old = x_;
  // Compute new coordinates.
  x_[0] = cost*x_old[0] + sint*x_old[2];
  // x_[1] is unchanged
  x_[2] = -sint*x_old[0] + cost*x_old[2];
}

void ThreeVector::rotate_around_z(double theta) {
  double cost = std::cos(theta);
  double sint = std::sin(theta);
  // Get original coordinates.
  std::array<double, 3> x_old = x_;
  // Compute new coordinates.
  x_[0] = cost*x_old[0] - sint*x_old[1];
  x_[1] = sint*x_old[0] + cost*x_old[1];
  // x_[2] is unchanged
}

void ThreeVector::rotate_to(ThreeVector &r) {
  rotate_around_y(r.get_theta());
  rotate_around_z(r.get_phi());
}

std::ostream &operator<<(std::ostream &out, const ThreeVector &v) {
  using namespace std;
  out.put('(');
  out.fill(' ');
  for (auto x : v) {
    out << field<8> << x;
  }
  return out << ')';
}

}  // namespace Smash
