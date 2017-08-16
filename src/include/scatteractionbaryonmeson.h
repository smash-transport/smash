/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SCATTERACTIONBARYONMESON_H_
#define SRC_INCLUDE_SCATTERACTIONBARYONMESON_H_

#include "scatteraction.h"

namespace Smash {

/**
 * \ingroup action
 * ScatterActionBaryonMeson is a special ScatterAction which represents the
 * scattering of a baryon and a meson.
 */
class ScatterActionBaryonMeson : public ScatterAction {
 public:
  /* Inherit constructor. */
  using ScatterAction::ScatterAction;
  /** Determine the parametrized total cross section at high energies
   * for a baryon-meson collision. */
  double high_energy_cross_section() const override;

 protected:
  /**
   * \ingroup logging
   * Writes information about this scatter action to the \p out stream.
   */
  void format_debug_output(std::ostream &out) const override;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_SCATTERACTIONBARYONMESON_H_
