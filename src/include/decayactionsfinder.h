/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_DECAYACTIONSFINDER_H_
#define SRC_INCLUDE_DECAYACTIONSFINDER_H_

#include "actionfinderfactory.h"

namespace Smash {

class DecayActionsFinder : public ActionFinderFactory {
 public:
  std::vector<ActionPtr> find_possible_actions(Particles *particles, const ExperimentParameters &parameters, CrossSections *cross_sections = NULL)
      const override;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_DECAYACTIONSFINDER_H_
