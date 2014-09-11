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

namespace Smash {

DecayActionsFinder::DecayActionsFinder(const ExperimentParameters &parameters)
                     : ActionFinderFactory(parameters.timestep_duration()) {
}

ActionList DecayActionsFinder::find_possible_actions(Particles *particles) const {
  ActionList actions;
  actions.reserve(10);  // for short time steps this seems reasonable to expect
                        // less than 10 decays in most time steps

  for (const auto &p : particles->data()) {
    if (p.type().is_stable()) {
      continue;      /* particle doesn't decay */
    }

    /* The clock goes slower in the rest frame of the resonance */
    double resonance_frame_timestep = dt_ * p.inverse_gamma();

    // Create a candidate Action
    DecayAction act(p);
    const auto width = act.weight();  // total decay width (mass-dependent)

    /* Exponential decay. Lifetime tau = 1 / width
     * t / tau = width * t (remember GeV-fm conversion)
     * P(decay at Delta_t) = width * Delta_t
     * P(alive after n steps) = (1 - width * Delta_t)^n
     * = (1 - width * Delta_t)^(t / Delta_t)
     * -> exp(-width * t) when Delta_t -> 0
     */
    if (Random::canonical() < resonance_frame_timestep * width / hbarc) {
      /* Time is up! Set the particle to decay at this timestep. */
      actions.emplace_back(new DecayAction(std::move(act)));
    }
  }
  return std::move(actions);
}

}  // namespace Smash
