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

#include <vector>

#include "clock.h"
#include "forwarddeclarations.h"

namespace Smash {

/**
 * \ingroup action
 * ActionFinderInterface is the abstract base class for all action finders,
 * i.e. objects which create action lists.
 */
class ActionFinderInterface {
 public:
  /**
   * Abstract function for finding actions, given a list of particles.
   *
   * \param search_list a list of particles where each pair needs to be tested
   *                    for possible interaction
   * \param dt duration of the current time step in fm/c
   * \return The function returns a list (std::vector) of Action objects that
   *         could possibly be executed in this time step.
   */
  virtual ActionList find_possible_actions(const ParticleList &search_list,
                                           float dt) const = 0;
  /**
   * Abstract function for finding actions, given two lists of particles,
   * a search list and a neighbors list.
   *
   * \param search_list a list of particles where each pair needs to be tested
   *                    for possible interaction
   * \param neighbors_list a list of particles that need to be tested against
   *                       particles in search_list for possible interaction
   * \param dt duration of the current time step in fm/c
   * \return The function returns a list (std::vector) of Action objects that
   *         could possibly be executed in this time step.
   */
  virtual ActionList find_possible_actions(const ParticleList &search_list,
                                           const ParticleList &neighbors_list,
                                           float dt) const = 0;

  /**
   * This abstract function finds 'final' actions
   * (for cleaning up at the end of the simulation).
   */
  virtual ActionList find_final_actions(const Particles &) const = 0;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_ACTIONFINDERFACTORY_H_
