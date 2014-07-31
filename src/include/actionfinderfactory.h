/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_ACTIONFINDERFACTORY_H_
#define SRC_INCLUDE_ACTIONFINDERFACTORY_H_

#include "forwarddeclarations.h"

namespace Smash {

/**
 * ActionFinderFactory is the abstract base class for all action finders,
 * i.e. objects which create action lists.
 */
class ActionFinderFactory {
 public:
  /** Initialize the finder with the given parameters. */
  ActionFinderFactory(float dt) : dt_(dt) {}
  /** Pure virtual function for finding actions, given a list of particles. */
  virtual ActionList find_possible_actions(Particles *particles) const = 0;

 protected:
  /** Timestep duration. */
  float dt_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_ACTIONFINDERFACTORY_H_
