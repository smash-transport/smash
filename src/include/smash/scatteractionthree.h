/*
 *
 *    Copyright (c) 2015-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SCATTERACTIONTHREE_H_
#define SRC_INCLUDE_SCATTERACTIONTHREE_H_

#include "action.h"

namespace smash {

class ScatterActionThree : public Action {
 public:
  ScatterActionThree(const ParticleData& in_part1, const ParticleData& in_part2,
                     const ParticleData& in_part3, const ParticleData& out_part,
                     double time);

  void generate_final_state() override;

  double get_total_weight() const override;

  double get_partial_weight() const override;

 protected:
  /*
   * \ingroup logging
   * Writes information about this action to the \p out stream.
   */
  void format_debug_output(std::ostream& out) const override;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SCATTERACTIONTHREE_H_
