/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_POTENTIALS_H_
#define SRC_INCLUDE_POTENTIALS_H_

#include <vector>

#include "configuration.h"
#include "forwarddeclarations.h"
#include "particledata.h"
#include "threevector.h"

namespace Smash {

/**
 * A class that stores parameters of potentials, calculates
 * potentials and their gradients. Potentials are responsible
 * for long-range interactions and stand in the left part of
 * Boltzmann equation. Short-range interactions are taken into
 * account in the right part of it - in the collision term.
 */
class Potentials {
 public:
  explicit Potentials(Configuration conf);
  ~Potentials();
  /** Evaluates potential at point r. Potential is always taken in the local
   * Eckart rest frame, but point r is in the computational frame.
   *
   * \param[in] r Arbitrary space point where potential is calculated
   * \param[in] plist List of all particles to be used in \f$j^{\mu}\f$
   *            calculation. If the distance between particle and calculation
   *            point r, \f$ |r-r_i| > 6 \sigma \f$ then particle input
   *            to density will be ignored.
   * \param[in] gs_sigma Width of the gaussian (\f$ \sigma \f$),
   *  which represents particle Wigner density.
   **/
  double potential(ThreeVector r, const ParticleList &plist, double gs_sigma);

  /** Evaluates potential gradient at point r. Potential is always taken in
   * the local Eckart rest frame, but point r is in the computational frame.
   *
   * \param[in] r Arbitrary space point where potential gradient is calculated
   * \param[in] plist List of all particles to be used in \f$j^{\mu}\f$
   *            calculation. If the distance between particle and calculation
   *            point r, \f$ |r-r_i| > 6 \sigma \f$ then particle input
   *            to density will be ignored.
   * \param[in] gs_sigma Width of the gaussian (\f$ \sigma \f$),
   *  which represents particle Wigner density.
   **/
  ThreeVector potential_gradient(ThreeVector r, const ParticleList &plist,
                                                             double gs_sigma);
 private:
  // Skyrme potential on/off
  bool use_skyrme_;

  // Parameters of skyrme potentials
  double skyrme_a_, skyrme_b_, skyrme_tau_;

  // Symmetry potential on/off
  bool use_symmetry_;

  // Parameters of symmetry potential
  double symmetry_s_;
};
}  // namespace Smash

#endif  // SRC_INCLUDE_POTENTIALS_H_
