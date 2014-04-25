/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/decayactionsfinder.h"

#include <algorithm>    // std::sort

namespace Smash {

std::vector<ActionPtr> DecayActionsFinder::find_possible_actions(Particles *particles, const ExperimentParameters &parameters, CrossSections *cross_sections) const
{
  std::vector<ActionPtr> actions;
  FourVector velocity_lrf;
  velocity_lrf.set_x0(1.0);

  for (const auto &p : particles->data()) {
    std::vector<int> in_part;
    int id = p.id();
    float width = particles->type(id).width();
    /* particle doesn't decay */
    if (width < 0.0)
      continue;
    /* local rest frame velocity */
    velocity_lrf.set_x1(p.momentum().x1() / p.momentum().x0());
    velocity_lrf.set_x2(p.momentum().x2() / p.momentum().x0());
    velocity_lrf.set_x3(p.momentum().x3() / p.momentum().x0());

    /* The clock goes slower in the rest frame of the resonance */
    double inverse_gamma = sqrt(velocity_lrf.Dot(velocity_lrf));
    double resonance_frame_timestep = parameters.eps * inverse_gamma;

    /* Exponential decay. Average lifetime t_avr = 1 / width
     * t / t_avr = width * t (remember GeV-fm conversion)
     * P(decay at Delta_t) = width * Delta_t
     * P(alive after n steps) = (1 - width * Delta_t)^n
     * = (1 - width * Delta_t)^(t / Delta_t)
     * -> exp(-width * t) when Delta_t -> 0
     */
    if (drand48() < resonance_frame_timestep * width / hbarc) {
      /* Time is up! Set the particle to decay at this timestep */
      in_part.push_back(id);
      actions.emplace_back(new DecayAction(in_part,0.,2));
    }
  }

  /* Sort action list by time. */
  std::sort (actions.begin(), actions.end());

  return actions;
}
  
}  // namespace Smash
