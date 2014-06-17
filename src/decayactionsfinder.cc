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
#include "include/crosssections.h"
#include "include/experimentparameters.h"
#include "include/fourvector.h"
#include "include/particles.h"
#include "include/random.h"

namespace Smash {

std::vector<ActionPtr> DecayActionsFinder::find_possible_actions(
    Particles *particles, const ExperimentParameters &parameters,
    CrossSections *) const {
  std::vector<ActionPtr> actions;
  FourVector velocity_lrf;

  for (const auto &p : particles->data()) {
    std::vector<int> in_part;
    int id = p.id();
    float width = p.total_width();
    /* particle doesn't decay */
    if (width <= 0.0) {
      continue;
    }
    /* local rest frame velocity */
    velocity_lrf = FourVector(1., p.velocity());

    /* The clock goes slower in the rest frame of the resonance */
    double inverse_gamma = sqrt(velocity_lrf.Dot(velocity_lrf));
    double resonance_frame_timestep = parameters.timestep_duration()
                                    * inverse_gamma;

    /* Exponential decay. Average lifetime t_avr = 1 / width
     * t / t_avr = width * t (remember GeV-fm conversion)
     * P(decay at Delta_t) = width * Delta_t
     * P(alive after n steps) = (1 - width * Delta_t)^n
     * = (1 - width * Delta_t)^(t / Delta_t)
     * -> exp(-width * t) when Delta_t -> 0
     */
    if (Random::canonical() < resonance_frame_timestep * width / hbarc) {
      /* Time is up! Set the particle to decay at this timestep */
      in_part.push_back(id);
      actions.emplace_back(new DecayAction(in_part,0.,2));
    }
  }
  return std::move(actions);
}

}  // namespace Smash
