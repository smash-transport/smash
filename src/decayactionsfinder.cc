/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/decayactionsfinder.h"

#include "smash/constants.h"
#include "smash/cxx14compat.h"
#include "smash/decayaction.h"
#include "smash/experimentparameters.h"
#include "smash/fourvector.h"
#include "smash/particles.h"
#include "smash/random.h"

namespace smash {

ActionList DecayActionsFinder::find_actions_in_cell(
    const ParticleList &search_list, double dt) const {
  ActionList actions;
  /* for short time steps this seems reasonable to expect
   * less than 10 decays in most time steps */
  actions.reserve(10);

  for (const auto &p : search_list) {
    if (p.type().is_stable()) {
      continue;  // particle doesn't decay
    }

    DecayBranchList processes =
        p.type().get_partial_widths_hadronic(p.effective_mass());
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
    const double decay_time = random::exponential<double>(
        /* The clock goes slower in the rest
         * frame of the resonance */
        one_over_hbarc * p.inverse_gamma() * width);
    /* If the particle is not yet formed at the decay time,
     * it should not be able to decay */
    if (decay_time < dt &&
        (p.formation_time() < (p.position().x0() + decay_time))) {
      /* => decay_time ∈ [0, dt[
       * => the particle decays in this timestep. */
      auto act = make_unique<DecayAction>(p, decay_time);
      act->add_decays(std::move(processes));
      actions.emplace_back(std::move(act));
    }
  }
  return actions;
}

ActionList DecayActionsFinder::find_final_actions(const Particles &search_list,
                                                  bool /*only_res*/) const {
  ActionList actions;

  for (const auto &p : search_list) {
    if (p.type().is_stable()) {
      continue;  // particle doesn't decay
    }
    auto act = make_unique<DecayAction>(p, 0.);
    act->add_decays(p.type().get_partial_widths(p.effective_mass()));
    actions.emplace_back(std::move(act));
  }
  return actions;
}

}  // namespace smash
