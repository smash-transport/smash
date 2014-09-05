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

#include "forwarddeclarations.h"

namespace Smash {

/**
 * \ingroup action
 * ActionFinderInterface is the abstract base class for all action finders,
 * i.e. objects which create action lists.
 */
class ActionFinderInterface {
 public:
  /** Initialize the finder with the given parameters. */
  ActionFinderInterface(float dt) : dt_(dt) {}

  /**
   * Abstract function for finding actions, given a list of particles.
   *
   * \param search_list a list of particles where each pair needs to be tested
   *                    for possible interaction
   * \param neighbors_list a list of particles that need to be tested against
   *                       particles in search_list for possible interaction
   * \param particles We should not need this parameter. But until the
   *                  interfaces are fixed, this is the pointer to the global
   *                  Particles map, which allows to modify particles globally.
   * \param parameters
   * \param cross_sections
   *
   * \return The function returns a list (std::vector) of Action objects that
   *         could possibly be executed in this time step.
   */
  virtual ActionList find_possible_actions(const ParticleList &search_list,
                                           const ParticleList &neighbors_list,
                                           const Particles &particles) const = 0;

 protected:
  /** Timestep duration. */
  float dt_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_ACTIONFINDERFACTORY_H_
