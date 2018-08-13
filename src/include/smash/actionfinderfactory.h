/*
 *
 *    Copyright (c) 2014-2018
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
#include "lattice.h"
#include "potentials.h"

namespace smash {

/**
 * \ingroup action
 * ActionFinderInterface is the abstract base class for all action finders,
 * i.e. objects which create action lists.
 */
class ActionFinderInterface {
 public:
  virtual ~ActionFinderInterface() = default;

  /**
   * Abstract function for finding actions, given a list of particles.
   *
   * \param[in] search_list a list of particles where each pair needs to be
   *                  tested  for possible interaction
   * \param[in] dt duration of the current time step [fm]
   * \return The function returns a list (std::vector) of Action objects that
   *         could possibly be executed in this time step.
   */
  virtual ActionList find_actions_in_cell(const ParticleList &search_list,
                                          double dt) const = 0;
  /**
   * Abstract function for finding actions, given two lists of particles,
   * a search list and a neighbors list.
   *
   * \param[in] search_list a list of particles where each particle needs to
   *                  be tested for possible interactions with the neighbors
   * \param[in] neighbors_list a list of particles that need to be tested
   *                  against particles in search_list for possible interaction
   * \param[in] dt duration of the current time step [fm]
   * \return The function returns a list (std::vector) of Action objects that
   *         could possibly be executed in this time step.
   */
  virtual ActionList find_actions_with_neighbors(
      const ParticleList &search_list, const ParticleList &neighbors_list,
      double dt) const = 0;

  /**
   * Abstract function for finding actions between a list of particles and
   * the surrounding particles.
   *
   * Note: The first list can be a subset of the second list.
   *
   * \param[in] search_list a list of particles where each particle needs to be
   *                    tested for possible interactions with the surrounding
   *                    particles
   * \param[in] surrounding_list a list of particles that need to be tested
   *                  against particles in search_list for possible interaction
   * \param[in] dt duration of the current time step [fm]
   * \return The function returns a list (std::vector) of Action objects that
   *         could possibly be executed in this time step.
   */
  virtual ActionList find_actions_with_surrounding_particles(
      const ParticleList &search_list, const Particles &surrounding_list,
      double dt) const = 0;

  /**
   * This abstract function finds 'final' actions (for cleaning up at the end
   * of the simulation, e.g. letting the remaining resonances decay).
   * \param[in] search_list a list of particles where each particle needs to be
   *                    tested for possible interactions with the surrounding
   *                    particles
   * \param[in] only_res this optional parameter requests that only actions
   *                 regarding resonances are considered (disregarding stable
   *                 particles)
   * \return The function returns a list (std::vector) of Action objects.
   */
  virtual ActionList find_final_actions(const Particles &search_list,
                                        bool only_res = false) const = 0;
};

}  // namespace smash

#endif  // SRC_INCLUDE_ACTIONFINDERFACTORY_H_
