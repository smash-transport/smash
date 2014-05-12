/*
 *    Copyright (c) 2012-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_FOURVECTOR_H_
#define SRC_INCLUDE_FOURVECTOR_H_

#include <array>
#include <cmath>

#include "include/threevector.h"

namespace Smash {

/**
 * The FourVector class holds relevant values in Minkowski spacetime
 * with (+, −, −, −) metric signature.
 *
 * The overloaded operators are build according to Andrew Koenig
 * recommendations where  the compound assignment operators is used as a
 * base for their non-compound counterparts. This means that the
 * operator + is implemented in terms of +=. The operator+ returns
 * a copy of it's result. + and friends are non-members, while
 * the compound assignment counterparts, changing the left
 * argument, are a member of the FourVector class.
 */
class FourVector {
 public:
  /// default constructor nulls the fourvector components
  FourVector() : x_{{0., 0., 0., 0.}} {}
  /// copy constructor
  FourVector(double y0, double y1, double y2, double y3)
      : x_{{y0, y1, y2, y3}} {}
  FourVector(double y0, ThreeVector vec)
      : x_{{y0, vec.x1(), vec.x2(), vec.x3()}} {}
  /* getters and setters for the four components */
  double inline x0(void) const;
  void inline set_x0(double t);
  double inline x1(void) const;
  void inline set_x1(double x);
  double inline x2(void) const;
  void inline set_x2(double y);
  double inline x3(void) const;
  void inline set_x3(double z);
  /// get the three-vector (spatial components)
  ThreeVector inline threevec() const;
  /// calculate the scalar product with another four-vector
  double inline Dot(const FourVector &a) const;
  /// calculate the square of the vector (which is a scalar)
  double inline sqr() const;
  /// calculate the absolute value
  double inline abs() const;
  /// calculate the square of the spatial three-vector
  double inline sqr3() const;
  /// calculate the absolute value of the spatial three-vector
  double inline abs3() const;
  /** Returns the FourVector boosted with velocity.
   *
   * The current FourVector is not changed.
   *
   * \param velocity (\f$\vec{v}\f$) is a ThreeVector representing the boost velocity
   */
  FourVector LorentzBoost(const ThreeVector &velocity) const;

  /* overloaded operators */
  bool inline operator==(const FourVector &a) const;
  bool inline operator!=(const FourVector &a) const;
  bool inline operator<(const FourVector &a) const;
  bool inline operator>(const FourVector &a) const;
  bool inline operator<=(const FourVector &a) const;
  bool inline operator>=(const FourVector &a) const;
  bool inline operator==(const double &a) const;
  bool inline operator!=(const double &a) const;
  bool inline operator<(const double &a) const;
  bool inline operator>(const double &a) const;
  bool inline operator<=(const double &a) const;
  bool inline operator>=(const double &a) const;
  FourVector inline operator+=(const FourVector &a);
  FourVector inline operator-=(const FourVector &a);
  FourVector inline operator*=(const double &a);
  FourVector inline operator/=(const double &a);

  using iterator = std::array<double, 4>::iterator;
  using const_iterator = std::array<double, 4>::const_iterator;

  /**
   * Returns an iterator starting at the 0th component.
   *
   * The iterator implements the RandomIterator concept. Thus, you can simply
   * write `begin() + 1` to get an iterator that points to the 1st component.
   */
  iterator begin() { return x_.begin(); }

  /**
   * Returns an iterator pointing after the 4th component.
   */
  iterator end() { return x_.end(); }

  /// const overload of the above
  const_iterator begin() const { return x_.begin(); }
  /// const overload of the above
  const_iterator end() const { return x_.end(); }

  /// \see begin
  const_iterator cbegin() const { return x_.cbegin(); }
  /// \see end
  const_iterator cend() const { return x_.cend(); }

 private:
  std::array<double, 4> x_;
};

double inline FourVector::x0(void) const {
  return x_[0];
}

void inline FourVector::set_x0(const double t) {
  x_[0] = t;
}

double inline FourVector::x1(void) const {
  return x_[1];
}

void inline FourVector::set_x1(const double x) {
  x_[1] = x;
}

double inline FourVector::x2(void) const {
  return x_[2];
}

void inline FourVector::set_x2(const double y) {
  x_[2] = y;
}

double inline FourVector::x3(void) const {
  return x_[3];
}

void inline FourVector::set_x3(const double z) {
  x_[3] = z;
}

ThreeVector inline FourVector::threevec() const {
  return ThreeVector(x_[1],x_[2],x_[3]);
}

/// check if all four vector components are equal
bool inline FourVector::operator==(const FourVector &a) const {
  return std::abs(x_[0] - a.x_[0]) < 1e-12 &&
         std::abs(x_[1] - a.x_[1]) < 1e-12 &&
         std::abs(x_[2] - a.x_[2]) < 1e-12 &&
         std::abs(x_[3] - a.x_[3]) < 1e-12;
}

/// use == operator for the inverse != check
bool inline FourVector::operator!=(const FourVector &a) const {
  return !(*this == a);
}

/// all four vector components are below comparison vector
bool inline FourVector::operator<(const FourVector &a) const {
  return (x_[0] < a.x_[0]) && (x_[1] < a.x_[1]) && (x_[2] < a.x_[2]) && (x_[3] < a.x_[3]);
}

/// use < operator for the inverse by switching arguments
bool inline FourVector::operator>(const FourVector &a) const {
  return a < *this;
}

/// use > operator for less equal
bool inline FourVector::operator<=(const FourVector &a) const {
  return !(*this > a);
}

/// use < operator for greater equal
bool inline FourVector::operator>=(const FourVector &a) const {
  return !(*this < a);
}

/// all vector components are equal to that number
bool inline FourVector::operator==(const double &a) const {
  return std::abs(x_[0] - a) < 1e-12 && std::abs(x_[1] - a) < 1e-12 &&
         std::abs(x_[2] - a) < 1e-12 && std::abs(x_[3] - a) < 1e-12;
}

/// use == operator for the inverse !=
bool inline FourVector::operator!=(const double &a) const {
  return !(*this == a);
}

/// all vector components are below that number
bool inline FourVector::operator<(const double &a) const {
  return (x_[0] < a) && (x_[1] < a) && (x_[2] < a) && (x_[3] < a);
}

/// all vector components are above that number
bool inline FourVector::operator>(const double &a) const {
  return (x_[0] > a) && (x_[1] > a) && (x_[2] > a) && (x_[3] > a);
}

/// all vector components are less equal that number
bool inline FourVector::operator<=(const double &a) const {
  return !(*this > a);
}

/// all vector components are greater equal that number
bool inline FourVector::operator>=(const double &a) const {
  return !(*this < a);
}

/// += assignement addition
FourVector inline FourVector::operator+=(const FourVector &a) {
  this->x_[0] += a.x_[0];
  this->x_[1] += a.x_[1];
  this->x_[2] += a.x_[2];
  this->x_[3] += a.x_[3];
  return *this;
}

/// addition +operator uses +=
inline FourVector operator+(FourVector a, const FourVector &b) {
  a += b;
  return a;
}

/// -= assignement subtraction
FourVector inline FourVector::operator-=(const FourVector &a) {
  this->x_[0] -= a.x_[0];
  this->x_[1] -= a.x_[1];
  this->x_[2] -= a.x_[2];
  this->x_[3] -= a.x_[3];
  return *this;
}

/// subtraction operator- uses -=
inline FourVector operator-(FourVector a, const FourVector &b) {
  a -= b;
  return a;
}

/// assignement factor multiplication
FourVector inline FourVector::operator*=(const double &a) {
  this->x_[0] *= a;
  this->x_[1] *= a;
  this->x_[2] *= a;
  this->x_[3] *= a;
  return *this;
}

/// factor multiplication uses *=
inline FourVector operator*(FourVector a, const double &b) {
  a *= b;
  return a;
}

/// assignement factor division
FourVector inline FourVector::operator/=(const double &a) {
  this->x_[0] /= a;
  this->x_[1] /= a;
  this->x_[2] /= a;
  this->x_[3] /= a;
  return *this;
}

/// factor division uses /=
inline FourVector operator/(FourVector a, const double &b) {
  a /= b;
  return a;
}

double inline FourVector::Dot(const FourVector &a) const {
  return x_[0] * a.x_[0] - x_[1] * a.x_[1] - x_[2] * a.x_[2] - x_[3] * a.x_[3];
}

double inline FourVector::sqr() const {
  return x_[0] * x_[0] - x_[1] * x_[1] - x_[2] * x_[2] - x_[3] * x_[3];
}

double inline FourVector::abs() const {
  return std::sqrt(this->sqr());
}

double inline FourVector::sqr3() const {
  return this->threevec().sqr();
}

double inline FourVector::abs3() const {
  return this->threevec().abs();
}

}  // namespace Smash

#endif  // SRC_INCLUDE_FOURVECTOR_H_
