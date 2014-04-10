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
#include "include/collisions.h"
#include "include/constants.h"
#include "include/experiment.h"
#include "include/outputroutines.h"

namespace Smash {

// check particle pairs for collision
void ModusDefault::check_collision_geometry(
    Particles *particles, CrossSections *cross_sections,
    std::list<int> *collision_list, size_t *rejection_conflict,
    const ExperimentParameters &parameters) {
  FourVector distance;
  double neighborhood_radius_squared =
      parameters.cross_section * fm2_mb * M_1_PI * 4;
  for (const ParticleData &data : particles->data()) {
    for (const ParticleData &data2 : particles->data()) {
      /* exclude check on same particle and double counting */
      if (data.id() >= data2.id()) {
        continue;
      }
      distance = data.position() - data2.position();
      /* skip particles that are double interaction radius length away
       * (3-product gives negative values
       * with the chosen sign convention for the metric)
       */
      if (-distance.DotThree() > neighborhood_radius_squared) {
        continue;
      }
      collision_criteria_geometry(particles, cross_sections, collision_list,
                                  parameters.eps, data.id(), data2.id(),
                                  rejection_conflict);
    }
  }
}

/*general propagation routine */
void ModusDefault::propagate(Particles *particles,
                             const ExperimentParameters &parameters) {
  FourVector distance, position;
  for (ParticleData &data : particles->data()) {
    /* propagation for this time step */
    distance.set_FourVector(parameters.eps,
                            data.velocity_x() * parameters.eps,
                            data.velocity_y() * parameters.eps,
                            data.velocity_z() * parameters.eps);
    printd("Particle %d motion: %g %g %g %g\n", data.id(), distance.x0(),
           distance.x1(), distance.x2(), distance.x3());
    position = data.position();
    position += distance;
    data.set_position(position);
  }
}

}  // namespace Smash
