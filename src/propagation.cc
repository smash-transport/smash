/*
 *
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 *
 */
#include "include/propagation.h"

#include <cstdio>
#include <vector>

#include "include/Box.h"
#include "include/Parameters.h"
#include "include/ParticleData.h"
#include "include/outputroutines.h"

/* propagate all particles */
void propagate_particles(std::vector<ParticleData> *particles,
  Parameters const &parameters, Box const &box) {
    FourVector distance, position;

    for (std::vector<ParticleData>::iterator i = particles->begin();
         i != particles->end(); ++i) {
      /* The particle has formed a resonance or has decayed
       * and is not active anymore
       */
      if (i->process_type() > 0)
        continue;

      /* propagation for this time step */
      distance.set_FourVector(parameters.eps(),
        i->velocity_x() * parameters.eps(),
        i->velocity_y() * parameters.eps(),
        i->velocity_z() * parameters.eps());
      printd("Particle %d motion: %g %g %g %g\n", i->id(),
         distance.x0(), distance.x1(), distance.x2(), distance.x3());

      /* treat the box boundaries */
      position = i->position();
      position += distance;
      position = boundary_condition(position, box);
      i->set_position(position);
      printd_position(*i);
    }
}
