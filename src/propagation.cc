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
#include <map>
#include <vector>

#include "include/Box.h"
#include "include/Parameters.h"
#include "include/Particles.h"
#include "include/ParticleData.h"
#include "include/ParticleType.h"
#include "include/outputroutines.h"

/* propagate all particles */
void propagate_particles(Particles *particles,
  Parameters const &parameters, Box const &box) {
    FourVector distance, position;

    for (std::map<int, ParticleData>::iterator i = particles->begin();
         i != particles->end(); ++i) {
      /* propagation for this time step */
      distance.set_FourVector(parameters.eps(),
        i->second.velocity_x() * parameters.eps(),
        i->second.velocity_y() * parameters.eps(),
        i->second.velocity_z() * parameters.eps());
      printd("Particle %d motion: %g %g %g %g\n", i->first,
         distance.x0(), distance.x1(), distance.x2(), distance.x3());

      /* treat the box boundaries */
      bool wall_hit = false;
      position = i->second.position();
      position += distance;
      position = boundary_condition(position, box, &wall_hit);
      if (wall_hit)
        write_oscar(particles->data(i->first), particles->type(i->first),
                    1, 1);
      i->second.set_position(position);
      if (wall_hit)
        write_oscar(particles->data(i->first), particles->type(i->first));
      printd_position(i->second);
    }
}
