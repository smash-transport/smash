/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_PAULIBLOCKING_H_
#define SRC_INCLUDE_PAULIBLOCKING_H_

#include "configuration.h"
#include "experimentparameters.h"
#include "forwarddeclarations.h"
#include "threevector.h"

namespace Smash {

/**
 * A class that stores parameters needed for Pauli blocking,
 * tabulates necessary integrals and computes phase-space
 * density. Pauli blocking is the way to go from classical
 * Boltzmann equation to quantum one, effectively reducing
 * cross-sections by \f$1 - f(r,p) \f$ factors, where
 * \f$f(r,p)\f$ is a phase-space density at coordinate
 * and momentum of a final-state fermion. Effective reduction
 * of cross-section is done via random rejection of reaction
 * with probability \f$ 1 - f\f$. More details can be found in
 * T. Gaitanos, A.B. Larionov, H. Lenske, U. Mosel,
 * Phys. Rev. C 81 (2010) 054316., section III B. Our implementation
 * mainly follows this article (and therefore GiBUU, see
 * gibuu.hepforge.org).
 */
class PauliBlocker {
 public:
  PauliBlocker(Configuration conf, const ExperimentParameters &parameters);
  ~PauliBlocker();

  // Returns phase-space density at the point (r,p)
  float phasespace_dens(const ThreeVector r, const ThreeVector p,
                        const Particles* particles) const;
 private:

  // Sigma of the gaussian used for smearing
  float sig_;

  // Radius, after which gaussians (used for averaging) are cut, fm
  float rc_;

  // Radius of averaging in coordinate space, fm
  float rr_;

  // Radius of averaging in momentum space, GeV/c
  float rp_;

  // Testparticles number
  int ntest_;
};
}  // namespace Smash

#endif  // SRC_INCLUDE_PAULIBLOCKING_H_
