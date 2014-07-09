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

#include <assert.h>
#include <sstream>

namespace Smash {

Action::Action(const ParticleList &in_part, float time_of_execution)
    : incoming_particles_(in_part), time_of_execution_(time_of_execution),
      total_weight_(0.) {}

Action::~Action() {}

float Action::weight() const {
  return total_weight_;
}

void Action::add_process (const ProcessBranch &p) {
  subprocesses_.push_back(p);
  total_weight_ += p.weight();
}

void Action::add_processes (const ProcessBranchList &pv) {
  for (const auto &proc : pv) {
    subprocesses_.push_back(proc);
    total_weight_ += proc.weight();
  }
}

bool Action::is_valid(const Particles &particles) const {
  for (const auto &part : incoming_particles_) {
    if (!particles.has_data(part.id())) {
      return false;
    }
  }
  return true;
}

ParticleList Action::incoming_particles() const {
  ParticleList l;
  for (const auto &part : incoming_particles_) {
    l.emplace_back(part);
  }
  return std::move(l);
}


ParticleList Action::choose_channel () {
  if (total_weight_ < really_small) {
    return ParticleList();
  }
  double random_interaction = Random::canonical();
  float interaction_probability = 0.0;
  /* Loop through all subprocesses and select one by Monte Carlo, based on
   * their weights.  */
  for (const auto &proc : subprocesses_) {
    if (proc.pdg_list().size() < 1
        || proc.pdg_list()[0] == PdgCode::invalid()) {
      continue;
    }
    interaction_probability += proc.weight() / total_weight_;
    if (random_interaction < interaction_probability) {
      return proc.particle_list();
    }
  }
  /* Should never get here. */
  std::stringstream ss;
  ss << "problem in choose_channel: " << subprocesses_.size() << " " <<
         interaction_probability << " " << total_weight_;
  throw std::runtime_error(ss.str());
}


void Action::check_conservation(const size_t &id_process) const {

  /* Check momentum conservation */
  FourVector momentum_difference;
  for (const auto &part : incoming_particles_) {
    momentum_difference += part.momentum();
  }
  for (const auto &p : outgoing_particles_) {
    momentum_difference -= p.momentum();
  }

  /* TODO: throw an exception */
  if (fabs(momentum_difference.x0()) > really_small) {
    printf("Process %zu\n", id_process);
    printf("Warning: E conservation violation %g\n",
           momentum_difference.x0());
  }
  if (fabs(momentum_difference.x1()) > really_small)
    printf("Warning: px conservation violation %g\n",
           momentum_difference.x1());
  if (fabs(momentum_difference.x2()) > really_small)
    printf("Warning: py conservation violation %g\n",
           momentum_difference.x2());
  if (fabs(momentum_difference.x3()) > really_small)
    printf("Warning: pz conservation violation %g\n",
           momentum_difference.x3());

  // TODO: check other conservation laws (baryon number etc)
}

}  // namespace Smash
