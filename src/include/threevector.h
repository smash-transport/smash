/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_THREEVECTOR_H_
#define SRC_INCLUDE_THREEVECTOR_H_

#include <array>

namespace Smash {

class ThreeVector {
 public:
  /// default constructor nulls the fourvector components
  ThreeVector() : x_{{0., 0., 0.}} {}
  /// copy constructor
  ThreeVector(double y1, double y2, double y3)
      : x_{{y1, y2, y3}} {}
  double inline x1(void) const;
  void inline set_x1(double x);
  double inline x2(void) const;
  void inline set_x2(double y);
  double inline x3(void) const;
  void inline set_x3(double z);
  double inline sqr(void) const;
  /* operators */
  ThreeVector inline operator- ();
  ThreeVector inline operator*= (const double &a);
  ThreeVector inline operator/= (const double &a);
 private:
  std::array<double, 3> x_;
};

double inline ThreeVector::x1(void) const {
  return x_[0];
}

void inline ThreeVector::set_x1(const double x) {
  x_[0] = x;
}

double inline ThreeVector::x2(void) const {
  return x_[1];
}

void inline ThreeVector::set_x2(const double y) {
  x_[1] = y;
}

double inline ThreeVector::x3(void) const {
  return x_[2];
}

void inline ThreeVector::set_x3(const double z) {
  x_[2] = z;
}

double inline ThreeVector::sqr(void) const {
  return x_[0] * x_[0] + x_[1] * x_[1] + x_[2] * x_[2];
}

ThreeVector inline ThreeVector::operator- () {
  x_[0] = -x_[0];
  x_[1] = -x_[1];
  x_[2] = -x_[2];
  return *this;
}

ThreeVector inline ThreeVector::operator*= (const double &a) {
  x_[0] *= a;
  x_[1] *= a;
  x_[2] *= a;
  return *this;
}

inline ThreeVector operator* (ThreeVector a, const double &b) {
  a *= b;
  return a;
}

ThreeVector inline ThreeVector::operator/= (const double &a) {
  x_[0] /= a;
  x_[1] /= a;
  x_[2] /= a;
  return *this;
}

ThreeVector inline operator/ (ThreeVector a, const double &b) {
  a /= b;
  return a;
}

}  // namespace Smash

#endif  // SRC_INCLUDE_THREEVECTOR_H_
