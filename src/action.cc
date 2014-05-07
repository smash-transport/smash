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

namespace Smash {

Action::Action(const std::vector<int> &in_part, float time_of_execution)
    : incoming_particles(in_part), time_of_execution_(time_of_execution),
      total_weight_(0.), interaction_type_(0) {}

Action::Action(const std::vector<int> &in_part, float time_of_execution,
               int interaction_type)
    : incoming_particles(in_part), time_of_execution_(time_of_execution),
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


}  // namespace Smash
