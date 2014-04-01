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

#include "particledata.h"

namespace Smash {

class Action {
 public:
  Action(const std::vector<int> &in_part, float time_of_execution, int interaction_type);
  Action(const std::vector<int> &in_part, float time_of_execution, int interaction_type, const std::vector<int> &out_part);
  virtual ~Action();

  /**
   * for sorting by time of execution
   */
  bool operator<(const Action &rhs) const {
    return time_of_execution_ < rhs.time_of_execution_;
  }

  /* Return the first and second incoming particle.  */
  int in1() const;
  int in2() const;

  int process_type(void) const;
  const std::vector<int> &final_state(void) const;

 private:
  std::vector<int> ingoing_particles_;
  float time_of_execution_;
  int interaction_type_;
  std::vector<int> outgoing_particles_;
};

using ActionPtr = std::unique_ptr<Action>;

}  // namespace Smash

#endif  // SRC_INCLUDE_ACTION_H_
