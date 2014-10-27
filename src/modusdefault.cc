/*
 *
 *    Copyright (c) 2013-2014
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
                             const OutputsList &) {
  const auto &log = logger<LogArea::ModusDefault>();
  FourVector distance, position;
  for (ParticleData &data : particles->data()) {
    /* propagation for this time step */
    distance = FourVector(0.0,
                          data.velocity() * parameters.timestep_duration());
    log.debug("Particle ", data, " motion: ", distance);
    position = data.position() + distance;
    position.set_x0(parameters.new_particle_time());
    data.set_4position(position);
  }
}



}  // namespace Smash
