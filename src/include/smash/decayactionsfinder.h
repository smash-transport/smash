/*
 *
 *    Copyright (c) 2014-2020,2022-2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_DECAYACTIONSFINDER_H_
#define SRC_INCLUDE_SMASH_DECAYACTIONSFINDER_H_

#include <vector>

#include "actionfinderfactory.h"
#include "input_keys.h"

namespace smash {

/**
 * \ingroup action
 * A simple decay finder:
 * Just loops through all particles and checks if they can decay during the next
 * timestep.
 */
class DecayActionsFinder : public ActionFinderInterface {
 public:
  /**
   * Initialize the finder
   *
   * \param[in] res_lifetime_factor The multiplicative factor to be applied to
   *                                resonance lifetimes; default is 1
   * \param[in] do_non_strong_decays whether to do non-strong decays at the end
   * \param[in] force_decays_at_end whether to enforce decays at the end
   */
  explicit DecayActionsFinder(double res_lifetime_factor,
                              bool do_non_strong_decays,
                              bool force_decays_at_end)
      : res_lifetime_factor_(res_lifetime_factor),
        do_final_non_strong_decays_(do_non_strong_decays),
        find_final_decays_(force_decays_at_end) {}

  /**
   * Check the whole particle list for decays.
   *
   * \param[in] search_list All particles in grid cell.
   * \param[in] dt Size of timestep [fm]
   * \return List with the found (Decay)Action objects.
   */
  ActionList find_actions_in_cell(
      const ParticleList &search_list, double dt, const double,
      const std::vector<FourVector> &) const override;

  /// Ignore the neighbor searches for decays
  ActionList find_actions_with_neighbors(
      const ParticleList &, const ParticleList &, double,
      const std::vector<FourVector> &) const override {
    return {};
  }

  /// Ignore the surrounding searches for decays
  ActionList find_actions_with_surrounding_particles(
      const ParticleList &, const Particles &, double,
      const std::vector<FourVector> &) const override {
    return {};
  }

  /**
   * Force all resonances to decay at the end of the simulation.
   *
   * \param[in] search_list All particles at the end of simulation.
   * \return List with the found (Decay)Action objects.
   */
  ActionList find_final_actions(const Particles &search_list) const override;

  /// Multiplicative factor to be applied to resonance lifetimes
  const double res_lifetime_factor_ = 1.;

  /// Do all non-strong decays (including weak and electro-magnetic ones)
  const bool do_final_non_strong_decays_;

  /// Whether to find final decay actions
  const bool find_final_decays_;

  /**
   * Whether to initial state particles can decay. Useful for analyzing
   * interactions involving one or more resonances.
   */
  const bool decay_initial_particles_ =
      InputKeys::collTerm_decayInitial.default_value();
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_DECAYACTIONSFINDER_H_
