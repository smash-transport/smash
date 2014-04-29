/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_ACTION_H_
#define SRC_INCLUDE_ACTION_H_

#include <vector>
#include <memory>

#include "particles.h"
#include "processbranch.h"

namespace Smash {

class Action {
 public:
  Action(const std::vector<int> &in_part, float time_of_execution);
  Action(const std::vector<int> &in_part, float time_of_execution, int interaction_type);
  virtual ~Action();

  /* For sorting by time of execution. */
  bool operator<(const Action &rhs) const {
    return time_of_execution_ < rhs.time_of_execution_;
  }

  /* Returns the total weight, which is a cross section in case of a ScatterAction
   * and a decay width in case of a DecayAction. */
  float weight(void) const;

  /* These functions add new subprocesses.  */
  void add_process (ProcessBranch p);
  void add_processes (std::vector<ProcessBranch> &pv);

  /* Actually perform the action, e.g. carry out a decay or scattering.  */
  virtual void perform (Particles *particles, size_t &id_process) = 0;

 protected:
  std::vector<int> ingoing_particles_;
  float time_of_execution_;
  std::vector<ProcessBranch> subprocesses_;  /* list of possible subprocesses  */
  float total_weight_;                       /* sum of all subprocess weights  */
  /* Type of interaction: 0=elastic collision, 1=resonance formation, 2=decay */
  int interaction_type_;
  std::vector<int> outgoing_particles_;      /* PDG codes of final-state particles  */
};


class DecayAction : public Action {
 public:
  DecayAction (const std::vector<int> &in_part, float time_of_execution,
               int interaction_type);
  /* Decide for a particular decay channel via Monte-Carlo
   * and set the outgoing_particles_ correspondingly.  */
  void decide (Particles *particles);
  void perform (Particles *particles, size_t &id_process);
 private:
  int resonance_decay (Particles *particles);
  int one_to_two (Particles *particles);
  int one_to_three (Particles *particles);
};


class ScatterAction : public Action {
 public:
  ScatterAction (const std::vector<int> &in_part, float time_of_execution);
  /* Decide for a particular subprocess via Monte-Carlo
   * and set the outgoing_particles_ correspondingly.  */
  void decide ();
  void perform (Particles *particles, size_t &id_process);
};


using ActionPtr = std::unique_ptr<Action>;
using ScatterActionPtr = std::unique_ptr<ScatterAction>;

}  // namespace Smash

#endif  // SRC_INCLUDE_ACTION_H_
