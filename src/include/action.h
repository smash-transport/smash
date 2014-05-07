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

/**
 * Action is the base class for a generic process that takes a number of
 * incoming particles and transforms them into any number of outgoing particles.
 * Currently such an action can be either a decay or a two-body collision
 * (see derived classes).
 */
class Action {
 public:
  /** Simple constructor. */
  Action(const std::vector<int> &in_part, float time_of_execution);
  /** Constructor with known interaction_type. */
  Action(const std::vector<int> &in_part, float time_of_execution, int interaction_type);
  /** Destructor. */
  virtual ~Action();

  /** For sorting by time of execution. */
  bool operator<(const Action &rhs) const {
    return time_of_execution_ < rhs.time_of_execution_;
  }

  /** Returns the total weight, which is a cross section in case of a ScatterAction
   * and a decay width in case of a DecayAction. */
  float weight(void) const;

  /** Add a new subprocess.  */
  void add_process (ProcessBranch p);
  /** Add several new subprocesses at once.  */
  void add_processes (std::vector<ProcessBranch> &pv);

  /** Actually perform the action, e.g. carry out a decay or scattering.  */
  virtual void perform (Particles *particles, size_t &id_process) = 0;

 protected:
  /** ID codes of incoming particles  */
  std::vector<int> incoming_particles_;
  /** time at which the action is supposed to be performed  */
  float time_of_execution_;
  /** list of possible subprocesses  */
  std::vector<ProcessBranch> subprocesses_;
  /** sum of all subprocess weights  */
  float total_weight_;
  /** Type of interaction: 0=elastic collision, 1=resonance formation, 2=decay */
  int interaction_type_;
  /** PDG codes of final-state particles  */
  std::vector<PdgCode> outgoing_particles_;
};

/**
 * DecayAction is a special action which takes one single particle in the
 * initial state and makes it decay into a number of daughter particles
 * (currently two or three).
 */
class DecayAction : public Action {
 public:
  /** Constructor. */
  DecayAction (const std::vector<int> &in_part, float time_of_execution,
               int interaction_type);
  /** Decide for a particular decay channel via Monte-Carlo
   * and set the outgoing_particles_ correspondingly.  */
  void choose_channel (Particles *particles);
  /** Carry out the action, i.e. do the decay. */
  void perform (Particles *particles, size_t &id_process);
 private:
  int resonance_decay (Particles *particles);
  int one_to_two (Particles *particles);
  int one_to_three (Particles *particles);
};

/**
 * ScatterAction is a special action which takes two incoming particles
 * and performs a scattering, producing one or more final-state particles.
 */
class ScatterAction : public Action {
 public:
  /** Constructor. */
  ScatterAction (const std::vector<int> &in_part, float time_of_execution);
  /** Decide for a particular final-state channel via Monte-Carlo
   * and set the outgoing_particles_ correspondingly.  */
  void choose_channel ();
  /** Carry out the action, i.e. do the scattering. */
  void perform (Particles *particles, size_t &id_process);
};


using ActionPtr = std::unique_ptr<Action>;
using ScatterActionPtr = std::unique_ptr<ScatterAction>;

inline std::vector<ActionPtr> &operator+=(std::vector<ActionPtr> &lhs,
                                          std::vector<ActionPtr> &&rhs) {
  if (lhs.size() == 0) {
    lhs = std::move(rhs);
  } else {
    lhs.insert(lhs.end(), std::make_move_iterator(rhs.begin()),
               std::make_move_iterator(rhs.end()));
  }
  return lhs;
}

}  // namespace Smash

#endif  // SRC_INCLUDE_ACTION_H_
