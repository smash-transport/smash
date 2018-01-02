/*
 *
 *    Copyright (c) 2015-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_DECAYACTION_H_
#define SRC_INCLUDE_DECAYACTION_H_

#include <utility>

#include "action.h"

namespace smash {

/**
 * \ingroup action
 * DecayAction is a special action which takes one single particle in the
 * initial state and makes it decay into a number of daughter particles
 * (currently two or three).
 */
class DecayAction : public Action {
 public:
  /**
   * Construct a DecayAction from a particle \p p.
   *
   * It does not initialize the list of possible decay processes. You need to
   * call add_processes after construction.
   *
   * \param[in] p The particle that should decay if the action is performed.
   * \param[in] time Time at which the action is supposed to take place
   */
  DecayAction(const ParticleData &p, double time);

  /** Add several new decays at once. */
  void add_decays(DecayBranchList pv);

  /** Add one new decay.*/
  void add_decay(DecayBranchPtr p);

  /** Generate the final state of the decay process.
   * Performs a decay of one particle to two or three particles.
   *
   * \throws InvalidDecay
   */
  void generate_final_state() override;

  std::pair<double, double> sample_masses() const override;

  double raw_weight_value() const override { return total_width_; }

  double partial_weight() const override { return partial_width_; }

  double total_width() const { return total_width_; }

  /**
   * \ingroup exception
   * Thrown when DecayAction is called to perform with 0 or more than 2
   * entries in outgoing_particles.
   */
  class InvalidDecay : public std::invalid_argument {
    using std::invalid_argument::invalid_argument;
  };

 protected:
  /**
   * \ingroup logging
   * Writes information about this decay action to the \p out stream.
   */
  void format_debug_output(std::ostream &out) const override;

  /** list of possible decays  */
  DecayBranchList decay_channels_;

  /** total decay width */
  double total_width_;

  /** partial decay width to the chosen outgoing channel */
  double partial_width_;

  /** angular momentum of the decay */
  int L_ = 0;

 protected:
  /**
   * Kinematics of a 1-to-3 decay process.
   *
   * Sample the masses and momenta of the decay products in the
   * center-of-momentum frame.
   */
  virtual void one_to_three();
};

}  // namespace smash

#endif  // SRC_INCLUDE_DECAYACTION_H_
