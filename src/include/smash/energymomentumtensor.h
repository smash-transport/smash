/*
 *
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_ENERGYMOMENTUMTENSOR_H_
#define SRC_INCLUDE_ENERGYMOMENTUMTENSOR_H_

#include <cmath>

#include <array>
#include <string>

#include "fourvector.h"
#include "particledata.h"

namespace smash {

/**
 * \ingroup data
 *
 * The EnergyMomentumTensor class represents a symmetric positive
 * semi-definite energy-momentum tensor \f$ T^{\mu \nu}\f$.
 *
 * \fpPrecision
 * Energy-momentum tensor is basically constructed from particle
 * momenta, which are instances of the \c FourVector class, their
 * components being \c double. Also boosting and going to Landau frame
 * require good precision, if the Landau frame velocity of the system
 * is close to the speed of light.
 */
class EnergyMomentumTensor {
 public:
  /// The energy-momentum tensor is symmetric and has 10 independent components
  typedef std::array<double, 10> tmn_type;
  /// Default constructor (nulls all components)
  EnergyMomentumTensor() { Tmn_.fill(0.); }

  /**
   * Copy constructor
   * \param[in] Tmn Components of the energy-momentum tensor
   */
  explicit EnergyMomentumTensor(const tmn_type &Tmn) {
    for (size_t i = 0; i < 10; i++) {
      Tmn_[i] = Tmn[i];
    }
  }

  /// Return ith component of the tensor.
  double operator[](std::size_t i) const { return Tmn_[i]; }

  /**
   * Access the index of component \f$ (\mu, \nu) \f$.
   * \param[in] mu \f$\mu\f$ is the row index (0 to 3)
   * \param[in] nu \f$\nu\f$ is the line index (0 to 3)
   */
  static std::int8_t tmn_index(std::int8_t mu, std::int8_t nu) {
    // clang-format off
    constexpr std::array<std::int8_t, 16> indices = {0, 1, 2, 3,
                                                     1, 4, 5, 6,
                                                     2, 5, 7, 8,
                                                     3, 6, 8, 9};
    // clang-format on
    if (mu < 4 && nu < 4 && mu >= 0 && nu >= 0) {
      return indices[mu + 4 * nu];
    } else {
      throw std::invalid_argument("Invalid indices: " + std::to_string(mu) +
                                  ", " + std::to_string(nu));
    }
  }

  /// Addition of \f$T^{\mu \nu}_0\f$ to tensor
  EnergyMomentumTensor inline operator+=(const EnergyMomentumTensor &Tmn0);
  /// Subtraction of \f$T^{\mu \nu}_0\f$ of tensor
  EnergyMomentumTensor inline operator-=(const EnergyMomentumTensor &Tmn0);
  /// Scaling of the tensor by scalar \f$a\f$
  EnergyMomentumTensor inline operator*=(double a);
  /// Division of the tensor by scalar \f$a\f$
  EnergyMomentumTensor inline operator/=(double a);

  /// Iterator over the components
  using iterator = tmn_type::iterator;
  /// Constant iterator over the components
  using const_iterator = tmn_type::const_iterator;

  /**
   * Find the Landau frame 4-velocity from energy-momentum tensor.
   * IMPORTANT: resulting 4-velocity is fourvector with LOWER index
   */
  FourVector landau_frame_4velocity() const;

  /**
   * Boost to a given 4-velocity.
   * IMPORTANT: boost 4-velocity is fourvector with LOWER index
   * \param[in] u 4-velocity vector
   */
  EnergyMomentumTensor boosted(const FourVector &u) const;

  /**
   * Given momentum of the particle adds \f$p^{\mu}p^{\mu}/p^0\f$
   * to the energy momentum tensor.
   * \param[in] mom Momentum 4-vector with upper index, as all 4-momenta in
   * SMASH
   */
  void add_particle(const FourVector &mom);
  /**
   * Same, but \f$ p^{\mu}p^{\mu}/p^0\f$ times factor is added.
   * \param[in] p Reference to the data information of the particle
   * \see ParticleData
   * \param[in] factor Usually a smearing factor
   */
  void add_particle(const ParticleData &p, double factor);

  /**
   * Returns an iterator starting at the (0,0) component.
   *
   * The iterator implements the randomIterator concept. Thus, you can simply
   * write `begin() + 1` to get an iterator that points to the 1st component.
   */
  iterator begin() { return Tmn_.begin(); }

  /**
   * Returns an iterator pointing after the (3,3)th component.
   */
  iterator end() { return Tmn_.end(); }

  /// Const overload of the above
  const_iterator begin() const { return Tmn_.begin(); }
  /// Const overload of the above
  const_iterator end() const { return Tmn_.end(); }

  /// \see begin
  const_iterator cbegin() const { return Tmn_.cbegin(); }
  /// \see end
  const_iterator cend() const { return Tmn_.cend(); }

 private:
  /**
   * The internal storage of the components.
   * Tensor has 16 components, but it is symmetric, so number
   * of independent components reduces to 10.
   */
  tmn_type Tmn_;
};

/**
 * \ingroup logging
 * Prints out 4x4 tensor to the output stream.
 * \param[in] out Location of output
 * \param[in] Tmn Energy-momentum tensor
 */
std::ostream &operator<<(std::ostream &, const EnergyMomentumTensor &);

EnergyMomentumTensor inline EnergyMomentumTensor::operator+=(
    const EnergyMomentumTensor &Tmn0) {
  for (size_t i = 0; i < 10; i++) {
    Tmn_[i] += Tmn0[i];
  }
  return *this;
}
/// Direct addition operator
EnergyMomentumTensor inline operator+(EnergyMomentumTensor a,
                                      const EnergyMomentumTensor &b) {
  a += b;
  return a;
}

EnergyMomentumTensor inline EnergyMomentumTensor::operator-=(
    const EnergyMomentumTensor &Tmn0) {
  for (size_t i = 0; i < 10; i++) {
    Tmn_[i] -= Tmn0[i];
  }
  return *this;
}
/// Direct subtraction operator
EnergyMomentumTensor inline operator-(EnergyMomentumTensor a,
                                      const EnergyMomentumTensor &b) {
  a -= b;
  return a;
}

EnergyMomentumTensor inline EnergyMomentumTensor::operator*=(const double a) {
  for (size_t i = 0; i < 10; i++) {
    Tmn_[i] *= a;
  }
  return *this;
}
/// Direct multiplication operator
inline EnergyMomentumTensor operator*(EnergyMomentumTensor a, const double b) {
  a *= b;
  return a;
}
/// Direct multiplication operator
inline EnergyMomentumTensor operator*(const double a, EnergyMomentumTensor b) {
  b *= a;
  return b;
}

EnergyMomentumTensor inline EnergyMomentumTensor::operator/=(const double a) {
  for (size_t i = 0; i < 10; i++) {
    Tmn_[i] /= a;
  }
  return *this;
}
/// Direct division operator
EnergyMomentumTensor inline operator/(EnergyMomentumTensor a, const double b) {
  a /= b;
  return a;
}

}  // namespace smash

#endif  // SRC_INCLUDE_ENERGYMOMENTUMTENSOR_H_
