/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/decayactionsfinder.h"


namespace Smash {

std::vector<ActionPtr> DecayActionsFinder::find_possible_actions(Particles *particles, const ExperimentParameters &parameters, CrossSections *cross_sections) const
{
  std::vector<ActionPtr> actions;
  FourVector velocity_lrf;
  velocity_lrf.set_x0(1.0);

  for (auto i = particles->begin(); i != particles->end(); ++i) {
    std::vector<int> in_part;
    int id = i->first;
    /* particle doesn't decay */
    if (particles->type(id).width() < 0.0)
      continue;
    /* local rest frame velocity */
    velocity_lrf.set_x1(i->second.momentum().x1() / i->second.momentum().x0());
    velocity_lrf.set_x2(i->second.momentum().x2() / i->second.momentum().x0());
    velocity_lrf.set_x3(i->second.momentum().x3() / i->second.momentum().x0());

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
    if (drand48() < resonance_frame_timestep
                    * particles->type(id).width() / hbarc) {
      /* Time is up! Set the particle to decay at this timestep */
      in_part.push_back(id);
      actions.push_back(ActionPtr(new Action(in_part,0.,2)));
    }
  }

  return actions;
}
  
}  // namespace Smash
