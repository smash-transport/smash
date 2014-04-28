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

class ActionFinderFactory {
 public:
  virtual void find_possible_actions (std::vector<ActionPtr> &actions,
				      Particles *particles,
				      const ExperimentParameters &parameters,
				      CrossSections *cross_sections = nullptr) const = 0;
 private:
};

}  // namespace Smash

#endif  // SRC_INCLUDE_ACTIONFINDERFACTORY_H_
