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
#include <assert.h>

namespace Smash {

Action::Action(const std::vector<int> &in_part, float time_of_execution)
    : incoming_particles_(in_part), time_of_execution_(time_of_execution),
      total_weight_(0.), interaction_type_(0) {}

Action::Action(const std::vector<int> &in_part, float time_of_execution,
               int interaction_type)
    : incoming_particles_(in_part), time_of_execution_(time_of_execution),
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

bool Action::is_valid(const Particles &particles) const {
  for (int id : incoming_particles_) {
    if (!particles.has_data(id)) {
      return false;
    }
  }
  return true;
}

ParticleList Action::incoming_particles(const Particles &particles) const {
  assert(is_valid(particles));
  ParticleList l;
  for (int id : incoming_particles_) {
    l.emplace_back(particles.data(id));
  }
  return std::move(l);
}

void Action::check_conservation(const Particles *particles,
                                const size_t &id_process) const {

  /* Check momentum conservation */
  FourVector momentum_difference;
  for (const auto &i : incoming_particles_) {
    momentum_difference += particles->data(i).momentum();
  }
  for (const auto &p : outgoing_particles_) {
    momentum_difference -= p.momentum();
  }

  /* TODO: throw an exception */
  if (fabs(momentum_difference.x0()) > really_small) {
    printf("Process %zu type %i\n", id_process, interaction_type_);
    printf("Warning: Interaction type %i E conservation violation %g\n",
           interaction_type_, momentum_difference.x0());
  }
  if (fabs(momentum_difference.x1()) > really_small)
    printf("Warning: Interaction type %i px conservation violation %g\n",
           interaction_type_, momentum_difference.x1());
  if (fabs(momentum_difference.x2()) > really_small)
    printf("Warning: Interaction type %i py conservation violation %g\n",
           interaction_type_, momentum_difference.x2());
  if (fabs(momentum_difference.x3()) > really_small)
    printf("Warning: Interaction type %i pz conservation violation %g\n",
           interaction_type_, momentum_difference.x3());

  // TODO: check other conservation laws (baryon number etc)
}

}  // namespace Smash
