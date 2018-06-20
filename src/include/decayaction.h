/*
 *
 *    Copyright (c) 2015-2018
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

  /**
   * Add several new decays at once.
   * \param[in] pv List of decays to be added.
   */
  void add_decays(DecayBranchList pv);

  /**
   * Add one new decay.
   * \param[in] p Decay to be added.
   */
  void add_decay(DecayBranchPtr p);

  /**
   * Generate the final state of the decay process.
   * Performs a decay of one particle to two or three particles.
   *
   * \throws InvalidDecay
   */
  void generate_final_state() override;

  /**
   * Sample the masses of the final particles
   * \returns Pair of sampled masses of particle 1 and 2
   */
  std::pair<double, double> sample_masses() const override;

  /// Return the total width of the decay process.
  double get_total_weight() const override { return total_width_; }

  /**
   * Get partial width of chosen channel
   * \return Partial width of chosen channel
   */
  double get_partial_weight() const override { return partial_width_; }

  /**
   * Get total decay width
   * \return Total width of decay
   */
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

  /// List of possible decays
  DecayBranchList decay_channels_;

  /// total decay width
  double total_width_;

  /// partial decay width to the chosen outgoing channel
  double partial_width_;

  /// Angular momentum of the decay
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
