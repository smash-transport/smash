/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_THREEVECTOR_H_
#define SRC_INCLUDE_THREEVECTOR_H_

#include <array>
#include <cmath>

#include "constants.h"

namespace smash {

/**
 * \ingroup data
 *
 * The ThreeVector class represents a physical three-vector
 * \f$ \vec{x} = (x_1,x_2,x_3)\f$
 * with the components \f$ x_1,x_2,x_3 \f$.
 * It is related to the classes FourVector and Angles,
 * both of which can be converted into a ThreeVector using
 * the 'threevec()' method.
 *
 * \fpPrecision
 * \li The \c ThreeVector class, like the \c FourVector class, uses \c double
 * for storage, calculations, and in the interface. This is necessary because
 * most particles are close to light speed and thus the low-order bits in the
 * mantissa of the momentum make large differences in energy. If they are
 * discarded by rounding to single-precision, e.g. boosting to/from the
 * center-of-mass frame breaks.
 * \li It might be sufficient for \c ThreeVector of other quantities to use
 * single-precision, though. This could be implemented by making \c ThreeVector
 * a
 * class template and use \c ThreeVector<double> for momenta and \c
 * ThreeVector<float> for the rest.
 */
class ThreeVector {
 public:
  /// default constructor (nulls all components)
  ThreeVector() : x_({0., 0., 0.}) {}

  /// copy constructor
  ThreeVector(double y1, double y2, double y3) : x_({y1, y2, y3}) {}

  /// access the component at offset \p i.
  double &operator[](std::size_t i) { return x_[i]; }
  /// const overload of the above.
  double operator[](std::size_t i) const { return x_[i]; }

  /// \return first component
  double inline x1() const;
  /// set first component
  void inline set_x1(double x);
  /// \return second component
  double inline x2() const;
  /// set second component
  void inline set_x2(double y);
  /// \return third component
  double inline x3() const;
  /// set third component
  void inline set_x3(double z);
  /// \return the square of the vector (which is a scalar)
  double inline sqr() const;
  /// \return the absolute value
  double inline abs() const;
  /// \return the azimuthal angle phi
  double inline get_phi() const;
  /// \return the polar angle theta
  double inline get_theta() const;
  /**
   * Rotate vector by the given Euler angles phi, theta, psi. If we
   * assume the standard basis x, y, z then this means applying the
   * matrix for a rotation of phi about z, followed by the matrix for
   * a rotation theta about the rotated x axis. Last, psi is a rotation
   * about the rotated z axis.
   * \param[in] phi angle by which the first rotation is done about z axis.
   * \param[in] theta angle by which the second rotation is done
   *        about the rotated x axis.
   * \param[in] psi angle by which the third rotation is done
   *        about the rotated z axis.
   *
   * Euler angles are used to make rotation of several (different) position
   * vectors belonging to one rigid body easy. A ThreeVector could be rotated
   * via only two angles, but then the angles for rotating a rigid body
   * consisting of multiple particles would require a different pair of rotation
   * angles for every position.
   */
  void inline rotate(double phi, double theta, double psi);
  /**
   * Rotate the vector around the y axis by the given angle theta.
   * \param[in] theta angle by which the rotation is done about y axis.
   */
  void inline rotate_around_y(double theta);
  /**
   * Rotate the vector around the z axis by the given angle theta.
   * \param[in] theta angle by which the rotation is done about z axis.
   */
  void inline rotate_around_z(double theta);
  /**
   * Rotate the z-axis onto the vector r.
   * \param[in] r direction in which new new z-axis is aligned
   */
  void inline rotate_z_axis_to(ThreeVector &r);
  /// negation: Returns \f$-\vec x\f$
  ThreeVector inline operator-() const;
  /// increase this vector by \f$\vec v: \vec x^\prime = \vec x + \vec v\f$
  ThreeVector inline operator+=(const ThreeVector &v);
  /// decrease this vector by \f$\vec v: \vec x^\prime = \vec x - \vec v\f$
  ThreeVector inline operator-=(const ThreeVector &v);
  /// scale this vector by \f$a: \vec x^\prime = a \cdot \vec x\f$
  ThreeVector inline operator*=(const double &a);
  /// divide this vector by \f$a: \vec x^\prime = \frac{1}{a} \cdot \vec x\f$
  ThreeVector inline operator/=(const double &a);

  /// \return whether the vector is identical to another vector
  bool operator==(const ThreeVector &rhs) const { return x_ == rhs.x_; }
  /// \return whether the vector is different from another vector
  bool operator!=(const ThreeVector &rhs) const { return x_ != rhs.x_; }

  /**
   * \return cross product of this vector and another vector
   * \f$ \vec{this} \times \vec{b} \f$
   */
  ThreeVector inline CrossProduct(const ThreeVector &b) const;
  /// iterates over the components
  using iterator = std::array<double, 3>::iterator;
  /// iterates over the components
  using const_iterator = std::array<double, 3>::const_iterator;

  /**
   * \return an iterator starting at the 0th component.
   *
   * The iterator implements the randomIterator concept. Thus, you can simply
   * write `begin() + 1` to get an iterator that points to the 1st component.
   */
  iterator begin() { return x_.begin(); }

  /// \return an iterator pointing after the 4th component.
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
  /// the internal storage of the components.
  std::array<double, 3> x_;
};

/**
 * \ingroup logging
 * Writes the three components of the vector to the output stream.
 */
std::ostream &operator<<(std::ostream &, const ThreeVector &);

double inline ThreeVector::x1() const { return x_[0]; }

void inline ThreeVector::set_x1(const double x) { x_[0] = x; }

double inline ThreeVector::x2() const { return x_[1]; }

void inline ThreeVector::set_x2(const double y) { x_[1] = y; }

double inline ThreeVector::x3() const { return x_[2]; }

void inline ThreeVector::set_x3(const double z) { x_[2] = z; }

ThreeVector inline ThreeVector::operator-() const {
  ThreeVector neg(-x_[0], -x_[1], -x_[2]);
  return neg;
}

ThreeVector inline ThreeVector::operator+=(const ThreeVector &v) {
  x_[0] += v.x_[0];
  x_[1] += v.x_[1];
  x_[2] += v.x_[2];
  return *this;
}

/// \return sum of two three-vectors \f$ \vec{a} + \vec{b} \f$.
ThreeVector inline operator+(ThreeVector a, const ThreeVector &b) {
  a += b;
  return a;
}

ThreeVector inline ThreeVector::operator-=(const ThreeVector &v) {
  x_[0] -= v.x_[0];
  x_[1] -= v.x_[1];
  x_[2] -= v.x_[2];
  return *this;
}

/// \return difference between two three-vectors \f$ \vec{a} - \vec{b} \f$.
ThreeVector inline operator-(ThreeVector a, const ThreeVector &b) {
  a -= b;
  return a;
}

ThreeVector inline ThreeVector::operator*=(const double &a) {
  x_[0] *= a;
  x_[1] *= a;
  x_[2] *= a;
  return *this;
}

/// multiply a three-vector by constant factor \f$ b \vec{a} \f$.
inline ThreeVector operator*(ThreeVector a, const double &b) {
  a *= b;
  return a;
}

/// multiply a three-vector by constant factor \f$ a \vec{b} \f$.
inline ThreeVector operator*(const double &a, ThreeVector b) {
  b *= a;
  return b;
}

/**
 * \return inner product of two three-vectors
 * \f$ \vec{a} \cdot \vec{b} \f$
 */
inline double operator*(ThreeVector a, const ThreeVector &b) {
  return a.x1() * b.x1() + a.x2() * b.x2() + a.x3() * b.x3();
}

/**
 * \return cross product of two three-vectors
 * \f$ \vec{a} \times \vec{b} \f$
 */
ThreeVector inline ThreeVector::CrossProduct(const ThreeVector &b) const {
  return ThreeVector(x_[1] * b.x3() - x_[2] * b.x2(),
                     x_[2] * b.x1() - x_[0] * b.x3(),
                     x_[0] * b.x2() - x_[1] * b.x1());
}

ThreeVector inline ThreeVector::operator/=(const double &a) {
  const double a_inv = 1.0 / a;
  x_[0] *= a_inv;
  x_[1] *= a_inv;
  x_[2] *= a_inv;
  return *this;
}

/// divide a three-vector by constant factor \f$ \vec{a} / b \f$.
ThreeVector inline operator/(ThreeVector a, const double &b) {
  a /= b;
  return a;
}

double inline ThreeVector::sqr() const { return (*this) * (*this); }

double inline ThreeVector::abs() const { return std::sqrt((*this) * (*this)); }

double inline ThreeVector::get_phi() const {
  if (std::abs(x1()) < really_small && std::abs(x2()) < really_small) {
    return 0.;
  } else {
    return std::atan2(x2(), x1());
  }
}

double inline ThreeVector::get_theta() const {
  double r = abs();
  return (r > 0.) ? std::acos(x3() / r) : 0.;
}

void inline ThreeVector::rotate(double phi, double theta, double psi) {
  // Compute the cosine and sine for each angle.
  const double cos_phi = std::cos(phi);
  const double sin_phi = std::sin(phi);
  const double cos_theta = std::cos(theta);
  const double sin_theta = std::sin(theta);
  const double cos_psi = std::cos(psi);
  const double sin_psi = std::sin(psi);
  // Get original coordinates.
  std::array<double, 3> x_old = x_;
  // Compute new coordinates.
  x_[0] = (cos_phi * cos_psi - sin_phi * cos_theta * sin_psi) * x_old[0] +
          (sin_phi * cos_psi + cos_phi * cos_theta * sin_psi) * x_old[1] +
          sin_theta * sin_psi * x_old[2];
  x_[1] = (-cos_phi * sin_psi - sin_phi * cos_theta * cos_psi) * x_old[0] +
          (-sin_phi * sin_psi + cos_phi * cos_theta * cos_psi) * x_old[1] +
          sin_theta * cos_psi * x_old[2];
  x_[2] = sin_phi * sin_theta * x_old[0] - cos_phi * sin_theta * x_old[1] +
          cos_theta * x_old[2];
}

void inline ThreeVector::rotate_around_y(double theta) {
  const double cost = std::cos(theta);
  const double sint = std::sin(theta);
  // Get original coordinates.
  std::array<double, 3> x_old = x_;
  // Compute new coordinates.
  x_[0] = cost * x_old[0] + sint * x_old[2];
  // x_[1] is unchanged
  x_[2] = -sint * x_old[0] + cost * x_old[2];
}

void inline ThreeVector::rotate_around_z(double theta) {
  const double cost = std::cos(theta);
  const double sint = std::sin(theta);
  // Get original coordinates.
  std::array<double, 3> x_old = x_;
  // Compute new coordinates.
  x_[0] = cost * x_old[0] - sint * x_old[1];
  x_[1] = sint * x_old[0] + cost * x_old[1];
  // x_[2] is unchanged
}

void inline ThreeVector::rotate_z_axis_to(ThreeVector &r) {
  rotate_around_y(r.get_theta());
  rotate_around_z(r.get_phi());
}

}  // namespace smash

#endif  // SRC_INCLUDE_THREEVECTOR_H_
