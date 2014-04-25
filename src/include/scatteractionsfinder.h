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
  /* Check for a single pair of particles (id_a, id_b) if a collision will happen
   * in the next timestep and create a corresponding Action object in that case. */
  ActionPtr check_collision(int id_a, int id_b, Particles *particles,
			    const ExperimentParameters &parameters,
			    CrossSections *cross_sections = nullptr) const;

  /* Check the whole particle list for collisions and return a list of Action objects. */
  std::vector<ActionPtr> find_possible_actions(Particles *particles,
					       const ExperimentParameters &parameters,
					       CrossSections *cross_sections = nullptr)
      const override;

 private:
};

}  // namespace Smash

#endif  // SRC_INCLUDE_SCATTERACTIONSFINDER_H_
