/*
 *
 *    Copyright (c) 2016
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SCATTERACTIONDELTAKAON_H_
#define SRC_INCLUDE_SCATTERACTIONDELTAKAON_H_

#include "scatteraction.h"
#include "scatteractionbaryonmeson.h"

namespace Smash {


/**
 * \ingroup action
 * ScatterActionDeltaKaon is a special ScatterActionBaryonMeson which represents the
 * scattering of a Delta and a kaon.
 */
class ScatterActionDeltaKaon : public ScatterActionBaryonMeson {
 public:
  /* Inherit constructor. */
  using ScatterActionBaryonMeson::ScatterActionBaryonMeson;
  /** Find all inelastic 2->2 processes for this reaction. */
  CollisionBranchList two_to_two_cross_sections() override;

 protected:
  /**
   * \ingroup logging
   * Writes information about this scatter action to the \p out stream.
   */
  void format_debug_output(std::ostream &out) const override;

 private:
  /**
   * Calculate cross sections for charge exchange in Delta-kaon
   * collisions. It is given by the cross section of the corresponding
   * Delta-kaon collision.
   */
  CollisionBranchList two_to_two_inel();
};

}  // namespace Smash

#endif  // SRC_INCLUDE_SCATTERACTIONDELTAKAON_H_
