/*
 *
 *    Copyright (c) 2014-2020,2022-2026
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
   * \param[in] par parameters from Experiment
   */
  explicit DecayActionsFinder(const ExperimentParameters &par)
      : res_lifetime_factor_(par.res_lifetime_factor),
        ignore_minimum_width_at_end_(par.ignore_minimum_width_at_end),
        force_decays_at_end_(par.force_decays_at_end),
        decay_initial_particles_(par.decay_initial_particles),
        spin_interaction_type_(par.spin_interaction_type) {}

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
  const double res_lifetime_factor_ =
      InputKeys::collTerm_resonanceLifetimeModifier.default_value();

  /// Do all non-strong decays (including weak and electro-magnetic ones)
  const bool ignore_minimum_width_at_end_ =
      InputKeys::collTerm_ignoreDecayWidthAtTheEnd.default_value();

  /// Whether to find final decay actions
  const bool force_decays_at_end_ =
      InputKeys::collTerm_forceDecaysAtEnd.default_value();

  /**
   * Whether to initial state particles can decay. Useful for analyzing
   * interactions involving one or more resonances.
   */
  const bool decay_initial_particles_ =
      InputKeys::collTerm_decayInitial.default_value();

  /// Spin interaction type
  const SpinInteractionType spin_interaction_type_ =
      InputKeys::collTerm_spinInteractions.default_value();
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_DECAYACTIONSFINDER_H_
