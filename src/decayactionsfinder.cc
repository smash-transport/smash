/*
 *
 *    Copyright (c) 2014-2020,2022-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/decayactionsfinder.h"

#include "smash/constants.h"
#include "smash/decayaction.h"
#include "smash/decaymodes.h"
#include "smash/fourvector.h"
#include "smash/random.h"

namespace smash {

ActionList DecayActionsFinder::find_actions_in_cell(
    const ParticleList &search_list, double dt, const double,
    const std::vector<FourVector> &) const {
  ActionList actions;
  /* for short time steps this seems reasonable to expect
   * less than 10 decays in most time steps */
  actions.reserve(10);

  for (const auto &p : search_list) {
    if (p.type().is_stable()) {
      continue;  // particle doesn't decay
    }

    if (!decay_initial_particles_ &&
        p.get_history().collisions_per_particle == 0) {
      continue;
    }

    DecayBranchList processes = p.type().get_partial_widths(
        p.momentum(), p.position().threevec(), WhichDecaymodes::Hadronic);
    // total decay width (mass-dependent)
    const double width = total_weight<DecayBranch>(processes);

    // check if there are any (hadronic) decays
    if (!(width > 0.0)) {
      continue;
    }

    constexpr double one_over_hbarc = 1. / hbarc;

    /* The decay_time is sampled from an exponential distribution.
     * Even though it may seem suspicious that it is sampled every
     * timestep, it can be proven that this still overall obeys
     * the exponential decay law.
     */
    double decay_time =
        res_lifetime_factor_ * random::exponential<double>(
                                   /* The clock goes slower in the rest
                                    * frame of the resonance */
                                   one_over_hbarc * p.inverse_gamma() * width);
    /* If the particle is not yet formed, shift the decay time by the time it
     * takes the particle to form */
    if (p.xsec_scaling_factor() < 1.0) {
      decay_time += p.formation_time() - p.position().x0();
    }
    if (decay_time < dt) {
      /* => decay_time ∈ [0, dt[
       * => the particle decays in this timestep. */
      auto act =
          std::make_unique<DecayAction>(p, decay_time, spin_interaction_type_);
      act->add_decays(std::move(processes));
      actions.emplace_back(std::move(act));
    }
  }
  return actions;
}

ActionList DecayActionsFinder::find_final_actions(
    const Particles &search_list) const {
  ActionList actions;
  if (find_final_decays_) {
    for (const auto &p : search_list) {
      if (!do_final_non_strong_decays_ && p.type().is_stable()) {
        continue;  // particle is stable with respect to strong interaction
      }

      if (p.type().decay_modes().is_empty()) {
        continue;  // particle cannot decay (not even e.m. or weakly)
      }

      auto act = std::make_unique<DecayAction>(p, 0., spin_interaction_type_);
      act->add_decays(p.type().get_partial_widths(
          p.momentum(), p.position().threevec(), WhichDecaymodes::All));
      actions.emplace_back(std::move(act));
    }
  }
  return actions;
}

}  // namespace smash
