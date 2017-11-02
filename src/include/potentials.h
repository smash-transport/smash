/*
 *
 *    Copyright (c) 2014-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_POTENTIALS_H_
#define SRC_INCLUDE_POTENTIALS_H_

#include <utility>
#include <vector>

#include "configuration.h"
#include "density.h"
#include "forwarddeclarations.h"
#include "particledata.h"
#include "threevector.h"

#ifdef BUILD_TESTS
#define VIRTUAL_FOR_TESTS virtual
#else
#define VIRTUAL_FOR_TESTS
#endif

namespace smash {

/**
 * A class that stores parameters of potentials, calculates
 * potentials and their gradients. Potentials are responsible
 * for long-range interactions and stand in the left part of
 * Boltzmann equation. Short-range interactions are taken into
 * account in the right part of it - in the collision term.
 */
class Potentials {
 public:
  Potentials(Configuration conf, const DensityParameters &parameters);
  ~Potentials();

  /// Evaluates skyrme potential given baryon density
  double skyrme_pot(const double baryon_density) const;
  /// Evaluates symmetry potential given baryon isospin density
  double symmetry_pot(const double baryon_isospin_density) const;

  /** Evaluates potential at point r. Potential is always taken in the local
   * Eckart rest frame, but point r is in the computational frame.
   *
   * \param[in] r Arbitrary space point where potential is calculated
   * \param[in] plist List of all particles to be used in \f$j^{\mu}\f$
   *            calculation. If the distance between particle and calculation
   *            point r, \f$ |r-r_i| > r_{cut} \f$ then particle input
   *            to density will be ignored.
   * \param[in] acts_on Type of particle on which potential is going to act
   *
   * \fpPrecision Why \c double?
   **/
  VIRTUAL_FOR_TESTS
  double potential(const ThreeVector &r, const ParticleList &plist,
                   const ParticleType &acts_on) const;

  /**
   * Evaluates the scaling factor of the forces acting on the particles. The
   * forces are equal to the product of the scaling factor and the gradient of
   * the potential. The scaling factors are usually less than one for hyperons,
   * and negative for anti-bayrons. The first component is the scaling factor
   * of the Skyrme force, and the second component is that of the symmetry
   * force.
   **/
  std::pair<double, int> force_scale(const ParticleType &data) const;

  /** Evaluates potential gradient at point r. Potential is always taken in
   * the local Eckart rest frame, but point r is in the computational frame.
   *
   * \param[in] r Arbitrary space point where potential gradient is calculated
   * \param[in] plist List of all particles to be used in \f$j^{\mu}\f$
   *            calculation. If the distance between particle and calculation
   *            point r, \f$ |r-r_i| > r_{cut} \f$ then particle input
   *            to density will be ignored.
   **/
  VIRTUAL_FOR_TESTS
  std::pair<ThreeVector, ThreeVector> potential_gradient(
      const ThreeVector &r, const ParticleList &plist) const;

  /// Is Skyrme potential on?
  VIRTUAL_FOR_TESTS
  bool use_skyrme() const { return use_skyrme_; }
  /// Is symmetry potential on?
  VIRTUAL_FOR_TESTS
  bool use_symmetry() const { return use_symmetry_; }

 private:
  /** Struct that contains gaussian sigma, cutoff and testparticle number
   *  needed for calculation
   */
  const DensityParameters param_;

  // Skyrme potential on/off
  bool use_skyrme_;

  // Symmetry potential on/off
  bool use_symmetry_;

  /** Parameters of skyrme potentials
   * \fpPrecision Why \c double?
   */
  double skyrme_a_, skyrme_b_, skyrme_tau_;

  /** Parameters of symmetry potential
   * \fpPrecision Why \c double?
   */
  double symmetry_s_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_POTENTIALS_H_
