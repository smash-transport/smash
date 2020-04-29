/*
 *
 *    Copyright (c) 2015-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SCATTERACTIONMULTI_H_
#define SRC_INCLUDE_SCATTERACTIONMULTI_H_

#include "action.h"

namespace smash {

class ScatterActionMulti : public Action {

  // TODO Write documentation

 public:
  ScatterActionMulti(const ParticleList& in_plist, double time);

  void generate_final_state() override;

  double get_total_weight() const override;

  double get_partial_weight() const override;

  void add_scattering();

  double probability_multi(double dt, const double cell_vol) const;

  /**
   * \ingroup exception
   * Thrown when ScatterAction is called to perform with unknown combination of
   * incoming and outgoing number of particles.
   */
  class InvalidScatterActionMulti : public std::invalid_argument {
    using std::invalid_argument::invalid_argument;
  };


 protected:
  /*
   * \ingroup logging
   * Writes information about this action to the \p out stream.
   */
  void format_debug_output(std::ostream& out) const override;

};

}  // namespace smash

#endif  // SRC_INCLUDE_SCATTERACTIONMULTI_H_
