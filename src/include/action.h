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

namespace Smash {

class Action {
 public:
  Action(const ParticleList &ingoing_particles_ids, float time_of_execution);
  virtual ~Action();

  /**
   * for sorting by time of execution
   */
  bool operator<(const Action &rhs) const {
    return time_of_execution_ < rhs.time_of_execution_;
  }

 private:
  ParticleList ingoing_particles_ids_;
  float time_of_execution_;
};

using ActionPtr = std::unique_ptr<Action>;

}  // namespace Smash

#endif  // SRC_INCLUDE_ACTION_H_
