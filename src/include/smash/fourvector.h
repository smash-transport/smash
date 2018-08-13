/*
 *    Copyright (c) 2012-2018
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

namespace smash {

/**
 * \ingroup data
 *
 * The FourVector class holds relevant values in Minkowski spacetime
 * with (+, −, −, −) metric signature.
 *
 * The overloaded operators are built according to Andrew Koenig
 * recommendations where the compound assignment operators is used as a
 * base for their non-compound counterparts. This means that the
 * operator + is implemented in terms of +=. The operator+ returns
 * a copy of its result. + and friends are non-members, while
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
  FourVector() : x_({0., 0., 0., 0.}) {}

  /**
   * copy constructor
   *
   * \param[in] y0 The time component to be copied
   * \param[in] y1 The x component to be copied
   * \param[in] y2 The y component to be copied
   * \param[in] y3 the z component to be copied
   */
  FourVector(double y0, double y1, double y2, double y3)
      : x_({y0, y1, y2, y3}) {}

  /**
   * construct from time-like component and a ThreeVector.
   *
   * \param[in] y0 The time component to be used
   * \param[in] vec A ThreeVector (x,y,z) to be used
   */
  FourVector(double y0, ThreeVector vec)
      : x_({y0, vec.x1(), vec.x2(), vec.x3()}) {}

  /**
   * access the component at offset \p i.
   * This operator results in the same as using the x0()...x3() functions
   *
   * \param[in] i the index of the component to access (has to be 0,1,2 or 3)
   * \return the component at index i
   */
  double &operator[](std::size_t i) { return x_[i]; }
  /// const overload of the [] operator
  double operator[](std::size_t i) const { return x_[i]; }

  /// \return the time-like component
  double inline x0() const;
  /// \param[in] t set time-like component
  void inline set_x0(double t);
  /// \return the first space-like component
  double inline x1() const;
  /// \param[in] x set first space-like component
  void inline set_x1(double x);
  /// \return the second space-like component
  double inline x2() const;
  /// \param[in] y set second space-like component
  void inline set_x2(double y);
  /// \return the third space-like component
  double inline x3() const;
  /// \param[in] z set third space-like component
  void inline set_x3(double z);
  /// \return the space-like three-vector (x,y,z components)
  ThreeVector inline threevec() const;

  /**
   * Get the velocity (3-vector divided by zero component).
   * Should only be used with momentum 4-vectors (not with space-time ones).
   *
   * \return the ThreeVector velocity
   */
  ThreeVector inline velocity() const;

  /**
   * calculate the scalar product with another four-vector
   *
   * \param[in] a the FourVector to dot product with *this
   * \return \f$x^\mu a_\mu\f$
   */
  double inline Dot(const FourVector &a) const;

  /**
   * calculate the square of the vector (which is a scalar)
   *
   * \return \f$x^\mu x_\mu\f$
   */
  double inline sqr() const;

  /**
   * calculate the lorentz invariant absolute value
   *
   * \return \f$\sqrt{x^\mu x_\mu}\f$
   *
   * Note that this will fail for space-like vectors.
   */
  double inline abs() const;

  /**
   * calculate the square of the spatial three-vector
   *
   * \return \f$\vec x \cdot \vec x\f$
   */
  double inline sqr3() const;

  /**
   * calculate the absolute value of the spatial three-vector
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
   *
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
   * x^\prime_i &=& x_i + v_i \cdot (\frac{\gamma - 1}{\vec{v}^2} \cdot
   * \vec{r}\cdot\vec{v} - \gamma \cdot x_0)\\
   *            &=& x_i + v_i \cdot (\frac{\gamma^2}{\gamma + 1} \cdot
   * \vec{r}\cdot\vec{v} - \gamma \cdot x_0)\\
   *            &=& x_i - v_i \cdot \frac{\gamma}{\gamma + 1} \cdot (\gamma(x_0
   * -
   * \vec{r}\cdot\vec{v}) + x_0 )\\
   *            &=& x_i - v_i \cdot \frac{\gamma}{\gamma + 1} \cdot (x^\prime_0
   * + x_0) \f}
   *
   * Note: This function is equivalent to -velocity Boost from ROOT
   */
  FourVector LorentzBoost(const ThreeVector &v) const;

  /**
   * Check if all four vector components are almost equal
   * (accuracy \f$10^{-4}\f$).
   *
   * \param[in] a The FourVector to compare to
   * \return Whether *this and a are almost equal
   */
  bool operator==(const FourVector &a) const;

  /**
   * checks inequality (logical complement to
   * FourVector::operator==(const FourVector&) const)
   *
   * \param[in] a The FourVector to compare to
   * \return Whether *this and a are not almost equal
   */
  bool inline operator!=(const FourVector &a) const;

  /**
   * checks if \f$x^\mu < a^\mu\f$ for all \f$\mu\f$
   * (all four vector components are below comparison vector)
   *
   * \param[in] a The FourVector to compare to
   * \return Whether all components of *this are strictly below the
   *         corresponding components of a
   */
  bool inline operator<(const FourVector &a) const;

  /**
   * checks if \f$x^\mu > a^\mu\f$ for all \f$\mu\f$
   * (all four vector components are above comparison vector)
   *
   * \param[in] a The FourVector to compare to
   * \return Whether all components of *this are strictly above the
   *         corresponding components of a
   */
  bool inline operator>(const FourVector &a) const;

  /**
   * logical complement to FourVector::operator>(const FourVector&) const
   *
   * \param[in] a The FourVector to compare to
   * \return Whether all components of *this are below or equal to the
   *         corresponding components of a
   */
  bool inline operator<=(const FourVector &a) const;

  /**
   * logical complement to FourVector::operator<(const FourVector&) const
   *
   * \param[in] a The FourVector to compare to
   * \return Whether all components of *this are above or equal to the
   *         corresponding components of a
   */
  bool inline operator>=(const FourVector &a) const;

  /**
   * adds \f$a_\mu: x_\mu^\prime = x_\mu + a_\mu\f$
   *
   * \param[in] a The FourVector to add
   * \return FourVector that consists of the added components of *this and a
   */
  FourVector inline operator+=(const FourVector &a);

  /**
   * subtracts \f$a_\mu: x_\mu^\prime = x_\mu - a_\mu\f$
   *
   * \param[in] a The FourVector to subtract
   * \return FourVector consisting of the components of a subtracted from *this
   */
  FourVector inline operator-=(const FourVector &a);

  /**
   * multiplies by \f$a: x_\mu^\prime = a \cdot x_\mu\f$
   *
   * \param[in] a The value with which to multiply
   * \return FourVector where each component of *this has been multiplied by a
   */
  FourVector inline operator*=(const double &a);

  /**
   * divides by \f$a: x_\mu^\prime = \frac{1}{a} \cdot x_\mu\f$
   *
   * \param[in] a The value by which to divide
   * \return FourVector where each component of *this has been divided by a
   */
  FourVector inline operator/=(const double &a);

  /// iterates over the components
  using iterator = std::array<double, 4>::iterator;
  /// iterates over the components
  using const_iterator = std::array<double, 4>::const_iterator;

  /**
   * \return An iterator starting at the 0th component.
   *
   * The iterator implements the randomIterator concept. Thus, you can simply
   * write `begin() + 1` to get an iterator that points to the 1st component.
   */
  iterator begin() { return x_.begin(); }

  /// \return An iterator pointing after the 4th component.
  iterator end() { return x_.end(); }

  /// \return A const_iterator pointing at the 0th component.
  const_iterator begin() const { return x_.begin(); }
  /// \return A const_iterator pointing after the 4th component.
  const_iterator end() const { return x_.end(); }

  /// \see begin
  const_iterator cbegin() const { return x_.cbegin(); }
  /// \see end
  const_iterator cend() const { return x_.cend(); }

 private:
  /// internal storage of this vector's components
  std::array<double, 4> x_;
};

// Definitions of previous inline functions

double inline FourVector::x0(void) const { return x_[0]; }

void inline FourVector::set_x0(const double t) { x_[0] = t; }

double inline FourVector::x1() const { return x_[1]; }

void inline FourVector::set_x1(const double x) { x_[1] = x; }

double inline FourVector::x2() const { return x_[2]; }

void inline FourVector::set_x2(const double y) { x_[2] = y; }

double inline FourVector::x3() const { return x_[3]; }

void inline FourVector::set_x3(const double z) { x_[3] = z; }

ThreeVector inline FourVector::threevec() const {
  return ThreeVector(x_[1], x_[2], x_[3]);
}

ThreeVector inline FourVector::velocity() const { return threevec() / x0(); }

// use == operator for the inverse != check
bool inline FourVector::operator!=(const FourVector &a) const {
  return !(*this == a);
}

bool inline FourVector::operator<(const FourVector &a) const {
  return (x_[0] < a.x_[0]) && (x_[1] < a.x_[1]) && (x_[2] < a.x_[2]) &&
         (x_[3] < a.x_[3]);
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

FourVector inline FourVector::operator+=(const FourVector &a) {
  this->x_[0] += a.x_[0];
  this->x_[1] += a.x_[1];
  this->x_[2] += a.x_[2];
  this->x_[3] += a.x_[3];
  return *this;
}

// addition +operator uses +=
/**
 * add two FourVectors
 *
 * \param[in] a The first FourVector to add
 * \param[in] b The second FourVector to add
 * \return \f$x^\mu = a^\mu + b^\mu\f$
 */
inline FourVector operator+(FourVector a, const FourVector &b) {
  a += b;
  return a;
}

FourVector inline FourVector::operator-=(const FourVector &a) {
  this->x_[0] -= a.x_[0];
  this->x_[1] -= a.x_[1];
  this->x_[2] -= a.x_[2];
  this->x_[3] -= a.x_[3];
  return *this;
}

// subtraction -operator uses -=
/**
 * subtract two FourVectors
 *
 * \param[in] a The FourVector from which to subtract
 * \param[in] b The FourVector to subtract
 * \return \f$x^\mu = a^\mu - b^\mu\f$
 */
inline FourVector operator-(FourVector a, const FourVector &b) {
  a -= b;
  return a;
}

FourVector inline FourVector::operator*=(const double &a) {
  this->x_[0] *= a;
  this->x_[1] *= a;
  this->x_[2] *= a;
  this->x_[3] *= a;
  return *this;
}

// factor multiplication uses *=
/**
 * multiply a vector with a scalar
 *
 * \param[in] a The FourVector to multiply
 * \param[in] b The value with which to multiply
 * \return \f$x^\mu = b \cdot a^\mu\f$
 */
inline FourVector operator*(FourVector a, double b) {
  a *= b;
  return a;
}
/**
 * multiply a vector with a scalar
 *
 * \param[in] b The value with which to multiply
 * \param[in] a The FourVector to multiply
 * \return \f$x^\mu = b \cdot a^\mu\f$
 */
inline FourVector operator*(double b, FourVector a) {
  a *= b;
  return a;
}

FourVector inline FourVector::operator/=(const double &a) {
  const double a_inv = 1.0 / a;
  this->x_[0] *= a_inv;
  this->x_[1] *= a_inv;
  this->x_[2] *= a_inv;
  this->x_[3] *= a_inv;
  return *this;
}

// factor division uses /=
/**
 * divide a vector by a scalar
 *
 * \param[in] a The FourVector to divide
 * \param[in] b The value with which to divide
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
  if (this->sqr() > -really_small) {
    return std::sqrt(std::abs(this->sqr()));
  } else {
    throw std::runtime_error("Absolute value of 4-vector could not be "
                    "determined, taking sqrt of negative value.");
  }
}

double inline FourVector::sqr3() const { return this->threevec().sqr(); }

double inline FourVector::abs3() const { return this->threevec().abs(); }

/**\ingroup logging
 * Writes the four components of the vector to the output stream.
 *
 * \param[in] os The ostream into which to output
 * \param[in] vec The FourVector to write into os
 */
std::ostream &operator<<(std::ostream &os, const FourVector &vec);

}  // namespace smash

#endif  // SRC_INCLUDE_FOURVECTOR_H_
