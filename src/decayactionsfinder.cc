/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/decayactionsfinder.h"

#include "include/constants.h"
#include "include/cxx14compat.h"
#include "include/decayaction.h"
#include "include/experimentparameters.h"
#include "include/fourvector.h"
#include "include/particles.h"
#include "include/random.h"

namespace Smash {

ActionList DecayActionsFinder::find_actions_in_cell(
    const ParticleList &search_list, float dt) const {
  ActionList actions;
  actions.reserve(10);  // for short time steps this seems reasonable to expect
                        // less than 10 decays in most time steps

  for (const auto &p : search_list) {
    if (p.type().is_stable()) {
      continue;      /* particle doesn't decay */
    }

    DecayBranchList processes =
                      p.type().get_partial_widths_hadronic(p.effective_mass());
    // total decay width (mass-dependent)
    const float width = total_weight<DecayBranch>(processes);

    // check if there are any (hadronic) decays
    if (!(width > 0.0)) {
      continue;
    }

    constexpr float one_over_hbarc = 1.f/static_cast<float>(hbarc);

    /* Exponential decay. Lifetime tau = 1 / width
     * t / tau = width * t (remember GeV-fm conversion)
     * P(decay at Delta_t) = width * Delta_t
     * P(alive after n steps) = (1 - width * Delta_t)^n
     * = (1 - width * Delta_t)^(t / Delta_t)
     * -> exp(-width * t) when Delta_t -> 0
     */
    const float decay_time = Random::exponential<float>(
        one_over_hbarc *
        p.inverse_gamma()  // The clock goes slower in the rest frame of the
                           // resonance
        * width);
    /* If the particle is not yet formed at the decay time,
     * it should not be able to decay */
    if (decay_time < dt && ((p.formation_time() > decay_time
        && p.cross_section_scaling_factor() > really_small) ||
        (p.formation_time() < decay_time))) {
      // => decay_time ∈ [0, dt[
      // => the particle decays in this timestep.
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
      continue;      /* particle doesn't decay */
    }
    auto act = make_unique<DecayAction>(p, 0.f);
    act->add_decays(p.type().get_partial_widths(p.effective_mass()));
    actions.emplace_back(std::move(act));
  }
  return actions;
}

}  // namespace Smash
