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

/**
 * \ingroup action
 * A simple decay finder:
 * Just loops through all particles and checks if they can decay during the next timestep.  */
class DecayActionsFinder : public ActionFinderFactory {
 public:
  /** Initialize the finder with the given parameters. */
  DecayActionsFinder(const ExperimentParameters &parameters);
  /** Check the whole particle list for decays
   * and return a list with the corrsponding Action objects. */
  ActionList find_possible_actions(Particles *particles) const override;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_DECAYACTIONSFINDER_H_
