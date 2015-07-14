/*
 *
 *    Copyright (c) 2013-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/propagation.h"

#include "include/logging.h"
#include "include/boxmodus.h"
#include "include/collidermodus.h"
#include "include/listmodus.h"
#include "include/spheremodus.h"

namespace Smash {

/* Simple straight line propagation without potentials*/
void propagate_straight_line(Particles *particles,
                             const ExperimentParameters &parameters) {
  const auto &log = logger<LogArea::Propagation>();
  const double dt = parameters.timestep_duration();
  for (ParticleData &data : *particles) {
    FourVector distance = FourVector(0.0, data.velocity() * dt);
    log.debug("Particle ", data, " motion: ", distance);
    FourVector position = data.position() + distance;
    position.set_x0(parameters.new_particle_time());
    data.set_4position(position);
    /* Test, if particle is formed and reset cross_section_scaling_factor
       TODO: Is there a way to only do that for particles that were unformed 
       in the previous timestep? */
    if (data.formation_time() < data.position().x0()) {
      data.set_cross_section_scaling_factor(1.0);
    }
  }
}

void propagate(Particles *particles, const ExperimentParameters &parameters,
               const Potentials &pot) {
  // Copy particles before propagation to calculate potentials from them
  const ParticleList plist = particles->copy_to_vector();
  const double dt = parameters.timestep_duration();
  const auto &log = logger<LogArea::Propagation>();

  for (ParticleData &data : *particles) {
    ThreeVector dU_dr = pot.potential_gradient(data.position().threevec(),
                                               plist, data.pdgcode());
    log.debug("Propagate: dU/dr = ", dU_dr);
    ThreeVector v = data.velocity();
    // predictor step assuming momentum-indep. potential, dU/dp = 0
    // then for momentum predictor = corrector
    data.set_4momentum(data.effective_mass(),
                       data.momentum().threevec() - dU_dr * dt);
    ThreeVector v_pred = data.velocity();
    // corrector step
    FourVector distance = FourVector(0.0, (v + v_pred) * (0.5 * dt));
    log.debug("Particle ", data, " motion: ", distance);
    FourVector position = data.position() + distance;
    position.set_x0(parameters.new_particle_time());
    data.set_4position(position);
  }
}

}  // namespace Smash
