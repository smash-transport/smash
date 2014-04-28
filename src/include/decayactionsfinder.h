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
  void find_possible_actions (std::vector<ActionPtr> &actions,
			      Particles *particles,
			      const ExperimentParameters &parameters,
			      CrossSections *cross_sections = nullptr) const override;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_DECAYACTIONSFINDER_H_
