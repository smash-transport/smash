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

#include "clock.h"
#include "forwarddeclarations.h"
#include <vector>

namespace Smash {

/**
 * \ingroup action
 * ActionFinderInterface is the abstract base class for all action finders,
 * i.e. objects which create action lists.
 */
class ActionFinderInterface {
 public:
  /** Initialize the finder with the given parameters. */
  ActionFinderInterface(Clock labclock) : labclock_(labclock) {}

  /**
   * Abstract function for finding actions, given a list of particles.
   *
   * \param search_list a list of particles where each pair needs to be tested
   *                    for possible interaction
   * \param neighbors_list a list of particles that need to be tested against
   *                       particles in search_list for possible interaction
   * \return The function returns a list (std::vector) of Action objects that
   *         could possibly be executed in this time step.
   */
  virtual ActionList find_possible_actions(
      const ParticleList &search_list,
      const std::vector<const ParticleList *> &neighbors_list) const = 0;

  /**
   * This abstract function finds 'final' actions
   * (for cleaning up at the end of the simulation).
   */
  virtual ActionList find_final_actions(const ParticleList &) const = 0;

 protected:
  /** Timestep duration. */
  Clock labclock_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_ACTIONFINDERFACTORY_H_
