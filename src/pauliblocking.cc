/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
*/

#include "include/pauliblocking.h"

namespace Smash {

/**
 * PauliBlocker constructor. Gets parameters from configuration.
 * Tabulates necessary integrals.
 */
PauliBlocker::PauliBlocker(Configuration conf, const ExperimentParameters &param)
    : sig_(param.gaussian_sigma),
      rc_(conf.take({"Rc"}, 2.2)),
      rr_(conf.take({"Rr"}, 1.86)),
      rp_(conf.take({"Rp"}, 0.08)),
      ntest_(param.testparticles) {

  /*!\Userguide
   * \page potentials PauliBlocker
   *
   * \key Rr (float, optional, default = 1.86): \n
   * Radius [fm] of sphere for averaging in the coordinate space
   *
   * \key Rc (float, optional, default = 2.2): \n
   * Radius [fm] at which gaussians used for smoothing are cut
   *
   * \key Rp (float, optional, default = 0.08): \n
   * Radius [GeV/c] of sphere for averaging in the momentum space
   */

}

PauliBlocker::~PauliBlocker() {
}

float PauliBlocker::phasespace_dens(const ThreeVector r, const ThreeVector p,
                                    const Particles* particles) const {
  // dummy return
  return 0.0;
};

}  // namespace Smash
