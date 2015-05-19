/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_ENERGYMOMENTUMTENSOR_H_
#define SRC_INCLUDE_ENERGYMOMENTUMTENSOR_H_

#include <array>
#include <cmath>

#include "fourvector.h"

namespace Smash {

/**
 * \ingroup data
 *
 * The EnergyMomentumTensor class represents a symmetric positive
 * semidifinite energy-momentum tensor \f$ T^{\mu \nu}\f$.
 *
 * \fpPrecision
 * Energy-momentum tensor is basically constructed from particle
 * momenta, which are instances of the \ref FourVector class, their
 * components being \c double. Also boosting and going to Landau frame
 * require good precision if the Landau frame velocity of the system
 * is close to the speed of light.
 */
class EnergyMomentumTensor {
 public:
  typedef std::array<double, 10> tmn_type;
  /// default constructor (nulls all components)
  EnergyMomentumTensor() {
    Tmn_.fill(0.);
  }

  /// copy constructor
  explicit EnergyMomentumTensor(const tmn_type& Tmn) {
    for (size_t i = 0; i < 10; i++) {
      Tmn_[i] = Tmn[i];
    }
  }

  const double &operator[](std::size_t i) { return Tmn_[i]; }
  double operator[](std::size_t i) const { return Tmn_[i]; }

  /// access the index of component \f$ (\mu, \nu) \f$.
  static std::int8_t tmn_index(std::int8_t mu, std::int8_t nu) {
    constexpr std::array<std::int8_t, 16> indices = {0, 1, 2, 3,
                                                     1, 4, 5, 6,
                                                     2, 5, 7, 8,
                                                     3, 6, 8, 9};
    return indices[mu + 4*nu];
    /* Another possibility (by lpang):
       inline inx(int i, int j){
         return (i<j) ? (7*i+2*j-i*i)>>1 : (7*j+2*i-j*j)>>1;
       }
    */
  }

  /// increase this tensor by \f$T^{\mu \nu}_0\f$
  EnergyMomentumTensor inline operator+= (const EnergyMomentumTensor &Tmn0);
  /// decrease this tensor by \f$T^{\mu \nu}_0f$
  EnergyMomentumTensor inline operator-= (const EnergyMomentumTensor &Tmn0);
  /// scale this tensor by scalar \f$a\f$
  EnergyMomentumTensor inline operator*= (double a);
  /// divide this tensor by scalar \f$a\f$
  EnergyMomentumTensor inline operator/= (double a);

  /// iterates over the components
  using iterator = tmn_type::iterator;
  /// iterates over the components
  using const_iterator = tmn_type::const_iterator;

  /**
   * Find the Landau frame 4-velocity from energy-momentum tensor.
   * IMPORTANT: resulting 4-velocity is fourvector with LOWER index
   */
  FourVector landau_frame_4velocity() const;

  /**
   * Boost to a given 4-velocity.
   * IMPORTANT: boost 4-velocity is fourvector with LOWER index
   */
  EnergyMomentumTensor boosted(const FourVector& u) const;

  /**
    * Given momentum p of the particle adds \f$ p^{\mu}p^{\mu}/p^0\f$
    * to the energy momentum tensor.
    * Input momentum is fourvector with upper index, as all 4-momenta in SMASH
    */
  void add_particle(const FourVector& mom);

  /**
   * Returns an iterator starting at the (0,0)th component.
   *
   * The iterator implements the RandomIterator concept. Thus, you can simply
   * write `begin() + 1` to get an iterator that points to the 1st component.
   */
  iterator begin() { return Tmn_.begin(); }

  /**
   * Returns an iterator pointing after the (3,3)th component.
   */
  iterator end() { return Tmn_.end(); }

  /// const overload of the above
  const_iterator begin() const { return Tmn_.begin(); }
  /// const overload of the above
  const_iterator end() const { return Tmn_.end(); }

  /// \see begin
  const_iterator cbegin() const { return Tmn_.cbegin(); }
  /// \see end
  const_iterator cend() const { return Tmn_.cend(); }

 private:
  /** The internal storage of the components.
    * Tensor has 16 components, but it is symmetric, so number
    * of independent components reduces to 10.
    */
  tmn_type Tmn_;
};

/**\ingroup logging
 * Prints out 4x4 tensor to the output stream.
 */
std::ostream &operator<<(std::ostream &, const EnergyMomentumTensor &);

EnergyMomentumTensor inline EnergyMomentumTensor::operator+= (
                                const EnergyMomentumTensor &Tmn0) {
  for (size_t i = 0; i < 10; i++) {
    Tmn_[i] += Tmn0[i];
  }
  return *this;
}

EnergyMomentumTensor inline operator+ (EnergyMomentumTensor a,
                                 const EnergyMomentumTensor &b) {
  a += b;
  return a;
}

EnergyMomentumTensor inline EnergyMomentumTensor::operator-= (
                                const EnergyMomentumTensor &Tmn0) {
  for (size_t i = 0; i < 10; i++) {
    Tmn_[i] -= Tmn0[i];
  }
  return *this;
}

EnergyMomentumTensor inline operator- (EnergyMomentumTensor a,
                                 const EnergyMomentumTensor &b) {
  a -= b;
  return a;
}

EnergyMomentumTensor inline EnergyMomentumTensor::operator*= (
                                               const double a) {
  for (size_t i = 0; i < 10; i++) {
    Tmn_[i] *= a;
  }
  return *this;
}

inline EnergyMomentumTensor operator* (EnergyMomentumTensor a,
                                               const double b) {
  a *= b;
  return a;
}

inline EnergyMomentumTensor operator* (const double a,
                                         EnergyMomentumTensor b) {
  b *= a;
  return b;
}

EnergyMomentumTensor inline EnergyMomentumTensor::operator/= (
                                                const double a) {
  for (size_t i = 0; i < 10; i++) {
    Tmn_[i] /= a;
  }
  return *this;
}

EnergyMomentumTensor inline operator/ (EnergyMomentumTensor a,
                                                const double b) {
  a /= b;
  return a;
}

}  // namespace Smash

#endif  // SRC_INCLUDE_ENERGYMOMENTUMTENSOR_H_
