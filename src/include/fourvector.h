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
#include <iosfwd>

#include "threevector.h"

namespace Smash {

/**
 * \ingroup data
 *
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
 *
 * \fpPrecision
 * \li The \c FourVector class uses \c double for storage, calculations, and in
 * the interface. This is necessary because most particles are close to light
 * speed and thus the low-order bits in the mantissa of the momentum make large
 * differences in energy. If they are discarded by rounding to single-precision,
 * e.g. boosting to/from the center-of-mass frame breaks. (see also
 * https://fias.uni-frankfurt.de/pm/projects/smash/wiki/Precision_considerations)
 * \n
 * \li It might be sufficient for \c FourVectors of other quantities to use
 * single-precision, though. This could be implemented by making \c FourVector a
 * class template and use \c FourVector<double> for momenta and \c
 * FourVector<float> for the rest.
 */
class FourVector {
 public:
  /// default constructor nulls the fourvector components
  FourVector() {
    x_ = {0., 0., 0., 0.};
  }
  /// copy constructor
  FourVector(double y0, double y1, double y2, double y3) {
    x_ = {y0, y1, y2, y3};
  }
  /// construct from time-like component and a ThreeVector.
  FourVector(double y0, ThreeVector vec) {
      x_ = {y0, vec.x1(), vec.x2(), vec.x3()};
  }

  /// access the component at offset \p i.
  double &operator[](std::size_t i) { return x_[i]; }
  /// const overload of the above.
  double operator[](std::size_t i) const { return x_[i]; }

  /// retrieve time-like component
  double inline x0() const;
  /// set time-like component
  void inline set_x0(double t);
  /// retrieve first space-like component
  double inline x1() const;
  /// set first space-like component
  void inline set_x1(double x);
  /// retrieve second space-like component
  double inline x2() const;
  /// set second space-like component
  void inline set_x2(double y);
  /// set third space-like component
  double inline x3() const;
  /// set third space-like component
  void inline set_x3(double z);
  /// get the three-vector (spatial components)
  ThreeVector inline threevec() const;
  /**
   * Get the velocity (3-vector divided by zero component).
   * Should only be used with momentum 4-vectors (not with space-time ones).
   */
  ThreeVector inline velocity() const;
  /** calculate the scalar product with another four-vector
   *
   * \return \f$x^\mu a_\mu\f$
   */
  double inline Dot(const FourVector &a) const;
  /** calculate the square of the vector (which is a scalar)
   *
   * \return \f$x^\mu x_\mu\f$
   */
  double inline sqr() const;
  /** calculate the lorentz invariant absolute value
   *
   * \return \f$\sqrt{x^\mu x_\mu}\f$
   *
   * Note that this will fail for space-like vectors.
   */
  double inline abs() const;
  /** calculate the square of the spatial three-vector
   *
   * \return \f$\vec x \cdot \vec x\f$
   *
   */
  double inline sqr3() const;
  /** calculate the absolute value of the spatial three-vector
   *
   * \return \f$\sqrt{\vec x \cdot \vec x}\f$
   */
  double inline abs3() const;
  /** Returns the FourVector boosted with velocity v.
   *
   * The current FourVector is not changed.
   *
   * \param v (\f$\vec{v}\f$) is a ThreeVector representing the
   * boost velocity
   */
  FourVector LorentzBoost(const ThreeVector &v) const;

  /// checks component-wise equality (accuracy \f$10^{-12}\f$)
  bool operator==(const FourVector &a) const;
  /// checks inequality (logical complement to
  /// FourVector::operator==(const FourVector&) const)
  bool inline operator!=(const FourVector &a) const;
  /// checks if \f$x^\mu > a^\mu\f$ for all \f$\mu\f$
  bool inline operator<(const FourVector &a) const;
  /// checks if \f$x^\mu > a^\mu\f$ for all \f$\mu\f$
  bool inline operator>(const FourVector &a) const;
  /// logical complement to FourVector::operator>(const FourVector&) const
  bool inline operator<=(const FourVector &a) const;
  /// logical complement to FourVector::operator<(const FourVector&) const
  bool inline operator>=(const FourVector &a) const;
  /// checks if \f$x^\mu < a\f$ for all \f$\mu\f$.
  bool inline operator==(const double &a) const;
  /// logical complement of FourVector::operator==(const double&) const
  bool inline operator!=(const double &a) const;
  /// checks if \f$x^\mu < a\f$ for all \f$\mu\f$.
  bool inline operator<(const double &a) const;
  /// checks if \f$x^\mu > a\f$ for all \f$\mu\f$.
  bool inline operator>(const double &a) const;
  /// logical complement to FourVector::operator>(const double&) const
  bool inline operator<=(const double &a) const;
  /// logical complement to FourVector::operator<(const double&) const
  bool inline operator>=(const double &a) const;
  /// adds \f$a_\mu: x_\mu^\prime = x_\mu + a_\mu\f$
  FourVector inline operator+=(const FourVector &a);
  /// subtracts \f$a_\mu: x_\mu^\prime = x_\mu - a_\mu\f$
  FourVector inline operator-=(const FourVector &a);
  /// multiplies by \f$a: x_\mu^\prime = a \cdot x_\mu\f$
  FourVector inline operator*=(const double &a);
  /// divides by \f$a: x_\mu^\prime = \frac{1}{a} \cdot x_\mu\f$
  FourVector inline operator/=(const double &a);

  /// iterates over the components
  using iterator = std::array<double, 4>::iterator;
  /// iterates over the components
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
  /// internal storage of this vector's components
  std::array<double, 4> x_;
};

double inline FourVector::x0(void) const {
  return x_[0];
}

void inline FourVector::set_x0(const double t) {
  x_[0] = t;
}

double inline FourVector::x1() const {
  return x_[1];
}

void inline FourVector::set_x1(const double x) {
  x_[1] = x;
}

double inline FourVector::x2() const {
  return x_[2];
}

void inline FourVector::set_x2(const double y) {
  x_[2] = y;
}

double inline FourVector::x3() const {
  return x_[3];
}

void inline FourVector::set_x3(const double z) {
  x_[3] = z;
}

ThreeVector inline FourVector::threevec() const {
  return ThreeVector(x_[1], x_[2], x_[3]);
}

ThreeVector inline FourVector::velocity() const {
  return threevec() / x0();
}

// use == operator for the inverse != check
bool inline FourVector::operator!=(const FourVector &a) const {
  return !(*this == a);
}

// all four vector components are below comparison vector
bool inline FourVector::operator<(const FourVector &a) const {
  return (x_[0] < a.x_[0]) && (x_[1] < a.x_[1]) &&
         (x_[2] < a.x_[2]) && (x_[3] < a.x_[3]);
}

// use < operator for the inverse by switching arguments
bool inline FourVector::operator>(const FourVector &a) const {
  return a < *this;
}

// use > operator for less equal
bool inline FourVector::operator<=(const FourVector &a) const {
  return !(*this > a);
}

// use < operator for greater equal
bool inline FourVector::operator>=(const FourVector &a) const {
  return !(*this < a);
}

// all vector components are equal to that number
bool inline FourVector::operator==(const double &a) const {
  return std::abs(x_[0] - a) < 1e-12 && std::abs(x_[1] - a) < 1e-12 &&
         std::abs(x_[2] - a) < 1e-12 && std::abs(x_[3] - a) < 1e-12;
}

// use == operator for the inverse !=
bool inline FourVector::operator!=(const double &a) const {
  return !(*this == a);
}

// all vector components are below that number
bool inline FourVector::operator<(const double &a) const {
  return (x_[0] < a) && (x_[1] < a) && (x_[2] < a) && (x_[3] < a);
}

// all vector components are above that number
bool inline FourVector::operator>(const double &a) const {
  return (x_[0] > a) && (x_[1] > a) && (x_[2] > a) && (x_[3] > a);
}

// all vector components are less equal that number
bool inline FourVector::operator<=(const double &a) const {
  return !(*this > a);
}

// all vector components are greater equal that number
bool inline FourVector::operator>=(const double &a) const {
  return !(*this < a);
}

// += assignement addition
FourVector inline FourVector::operator+=(const FourVector &a) {
  this->x_[0] += a.x_[0];
  this->x_[1] += a.x_[1];
  this->x_[2] += a.x_[2];
  this->x_[3] += a.x_[3];
  return *this;
}

// addition +operator uses +=
/** add two FourVectors
 *
 * \return \f$x^\mu = a^\mu + b^\mu\f$
 */
inline FourVector operator+(FourVector a, const FourVector &b) {
  a += b;
  return a;
}

// -= assignement subtraction
FourVector inline FourVector::operator-=(const FourVector &a) {
  this->x_[0] -= a.x_[0];
  this->x_[1] -= a.x_[1];
  this->x_[2] -= a.x_[2];
  this->x_[3] -= a.x_[3];
  return *this;
}

/** subtract two FourVectors
 *
 * \return \f$x^\mu = a^\mu - b^\mu\f$
 */
inline FourVector operator-(FourVector a, const FourVector &b) {
  a -= b;
  return a;
}

// assignement factor multiplication
FourVector inline FourVector::operator*=(const double &a) {
  this->x_[0] *= a;
  this->x_[1] *= a;
  this->x_[2] *= a;
  this->x_[3] *= a;
  return *this;
}

// factor multiplication uses *=
/** multiply a vector with a scalar
 *
 * \return \f$x^\mu = b \cdot a^\mu\f$
 */
inline FourVector operator*(FourVector a, double b) {
  a *= b;
  return a;
}
inline FourVector operator*(double b, FourVector a) {
  a *= b;
  return a;
}

// assignement factor division
FourVector inline FourVector::operator/=(const double &a) {
  this->x_[0] /= a;
  this->x_[1] /= a;
  this->x_[2] /= a;
  this->x_[3] /= a;
  return *this;
}

// factor division uses /=
/** divide a vector by a scalar
 *
 * \return \f$x^\mu = \frac{1}{b} \cdot a^\mu\f$
 */
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

/**\ingroup logging
 * Writes the four components of the vector to the output stream.
 */
std::ostream &operator<<(std::ostream &os, const FourVector &vec);

}  // namespace Smash

#endif  // SRC_INCLUDE_FOURVECTOR_H_
