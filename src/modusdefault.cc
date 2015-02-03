/*
 *
 *    Copyright (c) 2013-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include <cinttypes>
#include <list>

#include "include/modusdefault.h"
#include "include/experiment.h"
#include "include/logging.h"

namespace Smash {

// general propagation routine

void ModusDefault::propagate(Particles *particles,
                             const ExperimentParameters &parameters,
                             const OutputsList &, const Potentials* pot) {
  const auto &log = logger<LogArea::ModusDefault>();
  FourVector distance, position;
  const double dt = parameters.timestep_duration();
  if (!pot) {
    for (ParticleData &data : particles->data()) {
      /* propagation for this time step without potentials*/
      distance = FourVector(0.0, data.velocity() * dt);
      log.debug("Particle ", data, " motion: ", distance);
      position = data.position() + distance;
      position.set_x0(parameters.new_particle_time());
      data.set_4position(position);
    }
  } else {
    ThreeVector dU_dr, v, v_pred;
    ParticleList plist(particles->data().begin(), particles->data().end());

    for (ParticleData &data : particles->data()) {
      dU_dr = pot->potential_gradient(data.position().threevec(),
                                           plist, data.pdgcode());
      log.debug("Modusdef: dU/dr = ", dU_dr);
      v = data.velocity();
      // predictor step assuming momentum-indep. potential, dU/dp = 0
      // then for momentum predictor = corrector
      data.set_4momentum(data.effective_mass(),
                         data.momentum().threevec() - dU_dr * dt);
      v_pred = data.velocity();
      // corrector step
      distance = FourVector(0.0, (v + v_pred) * (0.5 * dt));

      log.debug("Particle ", data, " motion: ", distance);
      position = data.position() + distance;
      position.set_x0(parameters.new_particle_time());
      data.set_4position(position);
    }
  }
}

}  // namespace Smash
