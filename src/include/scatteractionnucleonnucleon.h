/*
 *
 *    Copyright (c) 2015-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SCATTERACTIONNUCLEONNUCLEON_H_
#define SRC_INCLUDE_SCATTERACTIONNUCLEONNUCLEON_H_

#include <utility>

#include "scatteractionbaryonbaryon.h"

namespace smash {

/**
 * \ingroup action
 * ScatterActionNucleonNucleon is a special ScatterAction which represents the
 * scattering of two nucleons.
 */
class ScatterActionNucleonNucleon : public ScatterActionBaryonBaryon {
 public:
  /* Inherit constructor. */
  using ScatterActionBaryonBaryon::ScatterActionBaryonBaryon;

 protected:
  /**
   * Sample final-state angles in a 2->2 collision (possibly anisotropic).
   *
   * \throws InvalidResonanceFormation
   */
  void sample_angles(std::pair<double, double> masses) override;

};

}  // namespace smash

#endif  // SRC_INCLUDE_SCATTERACTIONNUCLEONNUCLEON_H_
