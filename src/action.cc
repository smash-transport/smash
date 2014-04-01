/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/action.h"

namespace Smash {

Action::Action(const std::vector<int> &in_part, float time_of_execution,
	       int interaction_type)
    : ingoing_particles_(in_part), time_of_execution_(time_of_execution),
      interaction_type_(interaction_type) {}

Action::Action(const std::vector<int> &in_part, float time_of_execution,
	       int interaction_type, const std::vector<int> &out_part)
    : ingoing_particles_(in_part), time_of_execution_(time_of_execution),
      interaction_type_(interaction_type), outgoing_particles_(out_part) {}

Action::~Action() {}

// return first incoming particle
int Action::in1() const
{
  return ingoing_particles_[0];
}

// return second incoming particle
int Action::in2() const
{
  return ingoing_particles_[1];
}

/// return process type
int Action::process_type(void) const {
  return interaction_type_;
}

/// return final state
const std::vector<int> &Action::final_state(void) const {
  return outgoing_particles_;
}

}  // namespace Smash
