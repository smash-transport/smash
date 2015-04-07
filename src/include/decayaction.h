/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_DECAYACTION_H_
#define SRC_INCLUDE_DECAYACTION_H_

#include "action.h"

namespace Smash {


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
   * \param[in] time_of_execution time at which the action is supposed to take place
   */
  DecayAction(const ParticleData &p, float time_of_execution);

  /** Generate the final state of the decay process.
   * Performs a decay of one particle to two or three particles.
   *
   * \throws InvalidDecay
   */
  void generate_final_state() override;

  /**
   * \ingroup exception
   * Thrown when DecayAction is called to perform with 0 or more than 2
   * entries in outgoing_particles.
   */
  class InvalidDecay : public std::invalid_argument {
    using std::invalid_argument::invalid_argument;
  };

 protected:
  /// determine the total energy in the center-of-mass frame
  double sqrt_s() const override;

  /**
   * \ingroup logging
   * Writes information about this decay action to the \p out stream.
   */
  void format_debug_output(std::ostream &out) const override;

 private:
  /**
   * Kinematics of a 1-to-2 decay process.
   *
   * Sample the masses and momenta of the decay products in the
   * center-of-momentum frame.
   */
  void one_to_two();

  /**
   * Kinematics of a 1-to-3 decay process.
   *
   * Sample the masses and momenta of the decay products in the
   * center-of-momentum frame.
   */
  void one_to_three();
};


}  // namespace Smash

#endif  // SRC_INCLUDE_DECAYACTION_H_
