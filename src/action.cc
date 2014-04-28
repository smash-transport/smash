/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/action.h"

#include "include/constants.h"
#include "include/random.h"

namespace Smash {

Action::Action(const std::vector<int> &in_part, float time_of_execution)
    : ingoing_particles_(in_part), time_of_execution_(time_of_execution),
      total_weight_(0.), interaction_type_(0) {}

Action::Action(const std::vector<int> &in_part, float time_of_execution,
	       int interaction_type)
    : ingoing_particles_(in_part), time_of_execution_(time_of_execution),
      total_weight_(0.), interaction_type_(interaction_type) {}

Action::~Action() {}

float Action::weight(void) const {
  return total_weight_;
}

void Action::add_process (ProcessBranch p) {
  subprocesses_.push_back(p);
  total_weight_ += p.weight();
}

void Action::add_processes (std::vector<ProcessBranch> &pv) {
  for (auto proc = pv.begin(); proc != pv.end(); ++proc) {
    subprocesses_.push_back(*proc);
    total_weight_ += proc->weight();
  }
}

void Action::decide () {
  interaction_type_ = 0;
  if (total_weight_ > really_small) {
    double random_interaction = Random::canonical();
    float interaction_probability = 0.0;
    std::vector<ProcessBranch>::const_iterator proc = subprocesses_.begin();
    while (interaction_type_ == 0 && proc != subprocesses_.end()) {
      if (proc->particle_list().size() > 1
	  || proc->particle_list().at(0) != 0) {
	interaction_probability += proc->weight() / total_weight_;
	if (random_interaction < interaction_probability) {
	  interaction_type_ = proc->type();
	  outgoing_particles_ = proc->particle_list();
	}
      }
      ++proc;
    }
  }
}


}  // namespace Smash
