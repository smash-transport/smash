/*
 *
 *    Copyright (c) 2014-2017
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
#include "particles.h"
#include "pdgcode.h"
#include "threevector.h"

namespace smash {

/**
 * A class that stores parameters needed for Pauli blocking,
 * tabulates necessary integrals and computes phase-space
 * density. Pauli blocking is the way to go from the classical
 * Boltzmann equation to the quantum one, effectively reducing
 * cross-sections by \f$1 - f(r,p) \f$ factors, where
 * \f$f(r,p)\f$ is a phase-space density at coordinate
 * and momentum of a final-state fermion. Effective reduction
 * of cross-section is done via random rejection of reaction
 * with probability \f$ 1 - f\f$. More details can be found in
 * \iref{Gaitanos:2010fd}, section III B. Our implementation
 * mainly follows this article (and therefore GiBUU, see
 * http://gibuu.hepforge.org).
 */
class PauliBlocker {
 public:
  /**
   * PauliBlocker constructor. Gets parameters from configuration and
   * experiment. Tabulates necessary integrals.

   * \param[in] conf Configurations from config.yaml. 
   * \param[in] parameters Parameters given by Experiment.
   * \return The constructed object.
   */
  PauliBlocker(Configuration conf, const ExperimentParameters &parameters);

  /// Destructor
  ~PauliBlocker();

  /**
   * Calculate phase-space density of a particle species at the point (r,p).
   *
   * \param[in] r Position vector of the particle.
   * \param[in] p Momentum vector of the particle.
   * \param[in] particles List of all current particles.
   * \param[in] pdg PDG number of species for which density to be calculated.
   * \param[in] disregard Do not count particles that should be disregarded.
   *                       This is intended to avoid counting incoming 
   *                       particles when the phase-space density for outgoing
   *                       ones is estimated.
   * \return Phase-space density
   */
  double phasespace_dens(const ThreeVector &r, const ThreeVector &p,
                         const Particles &particles, const PdgCode pdg,
                         const ParticleList &disregard) const;

 private:
  /// Tabulate integrals for weights
  void init_weights();

  /// Analytical calculation of weights
  void init_weights_analytical();

  /// Standard deviation of the gaussian used for smearing
  double sig_;

  /// Radius, after which gaussians (used for averaging) are cut, fm
  double rc_;

  /// Radius of averaging in coordinate space, fm
  double rr_;

  /// Radius of averaging in momentum space, GeV/c
  double rp_;

  /// Testparticles number
  int ntest_;

  /// Weights: tabulated results of numerical integration
  std::array<double, 30> weights_;
};
}  // namespace smash

#endif  // SRC_INCLUDE_PAULIBLOCKING_H_
