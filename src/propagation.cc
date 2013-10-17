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

#include "include/Box.h"
#include "include/FourVector.h"
#include "include/Laboratory.h"
#include "include/Particles.h"
#include "include/ParticleData.h"
#include "include/outputroutines.h"

/* propagate all particles */
void propagate_particles(Particles *particles,
  Laboratory const &parameters, Box const &box) {
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

/* boundary_condition - enforce specific type of boundaries
 *
 * This assumes that the particle is at most one box length
 * away from the boundary to shift it in.
 */
FourVector boundary_condition(FourVector position, const Box &box,
                              bool *boundary_hit) {
  /* Check positivity and box size */
  if (position.x1() > 0 && position.x2() > 0 && position.x3() > 0
      && position.x1() < box.length() && position.x2() < box.length()
      && position.x3() < box.length())
    goto out;

  *boundary_hit = true;

  /* Enforce periodic boundary condition */
  if (position.x1() < 0)
    position.set_x1(position.x1() + box.length());

  if (position.x2() < 0)
    position.set_x2(position.x2() + box.length());

  if (position.x3() < 0)
    position.set_x3(position.x3() + box.length());

  if (position.x1() > box.length())
    position.set_x1(position.x1() - box.length());

  if (position.x2() > box.length())
    position.set_x2(position.x2() - box.length());

  if (position.x3() > box.length())
    position.set_x3(position.x3() - box.length());

 out:
    return position;
}
