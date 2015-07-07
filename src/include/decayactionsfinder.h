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

#include <vector>

#include "actionfinderfactory.h"

namespace Smash {

/**
 * \ingroup action
 * A simple decay finder:
 * Just loops through all particles and checks if they can decay during the next timestep.  */
class DecayActionsFinder : public ActionFinderInterface {
 public:
  /** Initialize the finder */
  DecayActionsFinder() {}
  /** Check the whole particle list for decays
   * and return a list with the corrsponding Action objects. */
  ActionList find_actions_in_cell(const ParticleList &search_list,
                                  float dt) const override;

  /// Ignore the neighbor searches for decays
  ActionList find_actions_with_neighbors(const ParticleList &,
                                         const ParticleList &,
                                         float) const override {
    return {};
  }
  /// Ignore the surrounding searches for decays
  ActionList find_actions_with_surrounding_particles(const ParticleList &,
                                                     const Particles &,
                                                     float) const override {
    return {};
  }

  /** Force all resonances to decay at the end of the simulation. */
  ActionList find_final_actions(const Particles &search_list) const override;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_DECAYACTIONSFINDER_H_
