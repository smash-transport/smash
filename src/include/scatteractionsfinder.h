/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SCATTERACTIONSFINDER_H_
#define SRC_INCLUDE_SCATTERACTIONSFINDER_H_

#include <vector>

#include "actionfinderfactory.h"

namespace Smash {

class ScatterActionsFinder : public ActionFinderFactory {
 public:
  std::vector<ActionPtr> find_possible_actions(Particles *particles,
					       const ExperimentParameters &parameters,
					       CrossSections *cross_sections = NULL)
      const override;

 private:
};

}  // namespace Smash

#endif  // SRC_INCLUDE_SCATTERACTIONSFINDER_H_
