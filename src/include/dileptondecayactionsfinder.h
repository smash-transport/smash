/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_DILEPTONDECAYACTIONSFINDER_H_
#define SRC_INCLUDE_DILEPTONDECAYACTIONSFINDER_H_

#include <vector>

#include "actionfinderfactory.h"

namespace Smash {

/**
 * \ingroup action
 * A simple dilepton decay finder:
 * Just loops through all particles and checks if they can decay into dileptons.
 */
class DileptonDecayActionsFinder : public ActionFinderInterface {
 public:
  /** Initialize the finder */
  DileptonDecayActionsFinder() {}
  /** Check the whole particle list for decays
   * and return a list with the corrsponding Action objects. */
  ActionList find_possible_actions(
      const ParticleList &search_list,
      const std::vector<const ParticleList *> &neighbors_list,
      float dt) const override;
  /** Force all resonances to decay at the end of the simulation. */
  ActionList find_final_actions(const Particles &search_list) const override;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_DILEPTONDECAYACTIONSFINDER_H_

