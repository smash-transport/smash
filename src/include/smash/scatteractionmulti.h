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
  // TODO Sensible (debug) logging for functions

 public:
  ScatterActionMulti(const ParticleList& in_plist, double time);

  void generate_final_state() override;

  double get_total_weight() const override;

  double get_partial_weight() const override;

  void add_possible_reactions(double dt, const double cell_vol);

  double probability() const { return total_probability_; }

  /**
   * \ingroup exception
   * TODO Thrown when ScatterAction is called to perform with unknown
   * combination of incoming and outgoing number of particles.
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

 private:
  void add_reaction(CollisionBranchPtr p);

  void add_reactions(CollisionBranchList pv);

  void annihilation();

  double probability_three_pi_to_one(const ParticleType& type_out, double dt,
                                       const double cell_vol) const;

  bool three_different_pions(const ParticleData& data_a,
                             const ParticleData& data_b,
                             const ParticleData& data_c) const;

  /// Total probability of reaction
  double total_probability_;

  /// Partial probability of the chosen outgoing channel
  double partial_probability_;

  /// List of possible collisions
  CollisionBranchList reaction_channels_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SCATTERACTIONMULTI_H_
