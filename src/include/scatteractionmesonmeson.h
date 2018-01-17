/*
 *
 *    Copyright (c) 2015-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SCATTERACTIONMESONMESON_H_
#define SRC_INCLUDE_SCATTERACTIONMESONMESON_H_

#include "scatteraction.h"

namespace smash {

/**
 * \ingroup action
 * ScatterActionMesonMeson is a special ScatterAction which represents the
 * scattering of two mesons.
 */
class ScatterActionMesonMeson : public ScatterAction {
 public:
  /* Inherit constructor. */
  using ScatterAction::ScatterAction;
  /**
   * Determine the (parametrized) hard non-diffractive string cross section
   * for a meson-meson collision.
   */
  double string_hard_cross_section() const override;

 protected:
  /**
   * \ingroup logging
   * Writes information about this scatter action to the \p out stream.
   */
  void format_debug_output(std::ostream &out) const override;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SCATTERACTIONMESONMESON_H_
