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

/** A simple decay finder:
 * Just loops through all particles and checks if they can decay during the next timestep.  */
class DecayActionsFinder : public ActionFinderFactory {
 public:
  /** Check the whole particle list for decays
   * and return a list with the corrsponding Action objects. */
  std::vector<ActionPtr> find_possible_actions(
      Particles *particles, const ExperimentParameters &parameters,
      CrossSections *cross_sections = nullptr) const override;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_DECAYACTIONSFINDER_H_
