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
  virtual std::vector<ActionPtr> find_possible_actions(Particles *particles,
						       CrossSections *cross_sections,
						       const ExperimentParameters &parameters)
      const = 0;

 private:
};

}  // namespace Smash

#endif  // SRC_INCLUDE_ACTIONFINDERFACTORY_H_
