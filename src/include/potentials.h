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
#include "threevector.h"
#include "particledata.h"
#include "forwarddeclarations.h"

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
  Potentials(Configuration conf);
  ~Potentials();

  double potential(ThreeVector r, const ParticleList &plist, double gs_sigma);
                            
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
