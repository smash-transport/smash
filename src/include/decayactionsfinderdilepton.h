/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_DECAYACTIONSFINDERDILEPTON_H_
#define SRC_INCLUDE_DECAYACTIONSFINDERDILEPTON_H_

#include <vector>

#include "actionfinderfactory.h"

namespace Smash {

/**
 * \ingroup action
 * A simple dilepton decay finder:
 * Just loops through all particles and checks if they can decay into dileptons.
 */
class DecayActionsFinderDilepton : public ActionFinderInterface {
 public:
  /** Initialize the finder */
  DecayActionsFinderDilepton() {}
  /** Check the whole particle list for decays
   * and return a list with the corrsponding Action objects. */
  ActionList find_possible_actions(
      const ParticleList &search_list,
      float dt) const override;

  /// Ignore the neighbor searches for (dilepton) decays
  ActionList find_possible_actions(const ParticleList &, const ParticleList &,
                                   float) const override {
    return {};
  }


  /** Force all resonances to decay at the end of the simulation. */
  ActionList find_final_actions(const Particles &search_list) const override;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_DECAYACTIONSFINDERDILEPTON_H_

