/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/decayactionsfinder.h"
#include "include/action.h"
#include "include/constants.h"
#include "include/experimentparameters.h"
#include "include/fourvector.h"
#include "include/particles.h"
#include "include/random.h"
#include "include/cxx14compat.h"

namespace Smash {

DecayActionsFinder::DecayActionsFinder(const ExperimentParameters &parameters)
    : ActionFinderInterface(parameters.timestep_duration()) {}

ActionList DecayActionsFinder::find_possible_actions(
    const ParticleList &search_list,
    const ParticleList &,  // the list of neighbors is irrelevant for decays
    const Particles &) const {
  ActionList actions;
  actions.reserve(10);  // for short time steps this seems reasonable to expect
                        // less than 10 decays in most time steps

  for (const auto &p : search_list) {
    if (p.type().is_stable()) {
      continue;      /* particle doesn't decay */
    }

    ProcessBranchList processes =
        p.type().get_partial_widths(p.effective_mass());
    const float width =
        total_weight(processes);  // total decay width (mass-dependent)

    /* Exponential decay. Lifetime tau = 1 / width
     * t / tau = width * t (remember GeV-fm conversion)
     * P(decay at Delta_t) = width * Delta_t
     * P(alive after n steps) = (1 - width * Delta_t)^n
     * = (1 - width * Delta_t)^(t / Delta_t)
     * -> exp(-width * t) when Delta_t -> 0
     *
     * A uniform distribution is not really correct, but good enough for small
     * time steps.
     */
    const float decay_time =
        Random::canonical<float>() *
        (static_cast<float>(hbarc) /
         (p.inverse_gamma()  // The clock goes slower in the rest frame of the
                             // resonance
          * width));

    if (decay_time < dt_) {
      // => decay_time ∈ [0, dt_[
      // => the particle decays in this timestep.
      auto act = make_unique<DecayAction>(p, decay_time);
      act->add_processes(std::move(processes));
      actions.emplace_back(std::move(act));
    }
  }
  return std::move(actions);
}

}  // namespace Smash
