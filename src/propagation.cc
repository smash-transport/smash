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

#include <vector>

#include "include/Box.h"
#include "include/Parameters.h"
#include "include/ParticleData.h"
#include "include/outputroutines.h"

/* propagate all particles */
void propagate_particles(std::vector<ParticleData> *particles,
  Parameters const &parameters, Box const &box) {
    FourVector distance, position;

    for (size_t i = 0; i < particles->size(); i++) {
      distance.set_FourVector(1.0, (*particles)[i].velocity_x(),
        (*particles)[i].velocity_y(), (*particles)[i].velocity_z());
      distance *= parameters.eps();
      printd("Particle %d motion: %g %g %g %g\n", (*particles)[i].id(),
         distance.x0(), distance.x1(), distance.x2(), distance.x3());

      /* treat the box boundaries */
      position = (*particles)[i].position();
      position += distance;
      position = boundary_condition(position, box);
      (*particles)[i].set_position(position);
      printd_position((*particles)[i]);
    }
}
