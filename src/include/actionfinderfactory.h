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

#include "action.h"
#include "particles.h"
#include "crosssections.h"
#include "experimentparameters.h"

namespace Smash {

/**
 * ActionFinderFactory is the abstract base class for all action finders,
 * i.e. objects which create action lists.
 */
class ActionFinderFactory {
 public:
  /** Pure virtual function for finding actions, given a list of particles. */
  virtual std::vector<ActionPtr> find_possible_actions(
      Particles *particles, const ExperimentParameters &parameters,
      CrossSections *cross_sections = nullptr) const = 0;

 private:
};

}  // namespace Smash

#endif  // SRC_INCLUDE_ACTIONFINDERFACTORY_H_
