/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SCATTERACTIONHYPERONPION_H_
#define SRC_INCLUDE_SCATTERACTIONHYPERONPION_H_

#include "scatteraction.h"
#include "scatteractionbaryonmeson.h"

namespace Smash {

/**
 * \ingroup action
 * ScatterActionHyperonPion is a special ScatterActionBaryonMeson which
 * represents the
 * scattering of a hyperon and a pion.
 */
class ScatterActionHyperonPion : public ScatterActionBaryonMeson {
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
};

}  // namespace Smash

#endif  // SRC_INCLUDE_SCATTERACTIONHYPERONPION_H_
