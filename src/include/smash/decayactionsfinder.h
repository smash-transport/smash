/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_DECAYACTIONSFINDER_H_
#define SRC_INCLUDE_DECAYACTIONSFINDER_H_

#include <vector>

#include "actionfinderfactory.h"

namespace smash {

/**
 * \ingroup action
 * A simple decay finder:
 * Just loops through all particles and checks if they can decay during the next
 * timestep.
 */
class DecayActionsFinder : public ActionFinderInterface {
 public:
  /// Initialize the finder
  DecayActionsFinder() {}

  /**
   * Check the whole particle list for decays.
   *
   * \param[in] search_list All particles in grid cell.
   * \param[in] dt Size of timestep [fm]
   * \return List with the found (Decay)Action objects.
   */
  ActionList find_actions_in_cell(const ParticleList &search_list,
                                  double dt) const override;

  /// Ignore the neighbor searches for decays
  ActionList find_actions_with_neighbors(const ParticleList &,
                                         const ParticleList &,
                                         double) const override {
    return {};
  }

  /// Ignore the surrounding searches for decays
  ActionList find_actions_with_surrounding_particles(const ParticleList &,
                                                     const Particles &,
                                                     double) const override {
    return {};
  }

  /**
   * Force all resonances to decay at the end of the simulation.
   *
   * \param[in] search_list All particles at the end of simulation.
   * \param[in] only_res optional parameter that requests that only actions
   *                     regarding resonances are considered (disregarding
   *                     stable particles) 
   * \return List with the found (Decay)Action objects.
   */
  ActionList find_final_actions(const Particles &search_list,
                                bool only_res = false) const override;
};

}  // namespace smash

#endif  // SRC_INCLUDE_DECAYACTIONSFINDER_H_
