/*
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#include <cstdio>

#include "../include/Particles.h"
#include "../include/constants.h"
#include "../include/ParticleData.h"
#include "../include/outputroutines.h"

int main() {
  /* checks for geometric distance criteria */
  ParticleData particle_a, particle_b;
  particle_a.set_id(0);
  particle_b.set_id(1);

  /* 2 particles with null momenta */
  particle_a.set_momentum(0.1, 0.0, 0.0, 0.0);
  particle_b.set_momentum(0.1, 0.0, 0.0, 0.0);
  particle_a.set_position(1.0, 1.0, 1.0, 1.0);
  particle_b.set_position(2.0, 2.0, 2.0, 2.0);

  /* check return of particle distance of null momenta particles */
  double distance_squared = particle_distance(&particle_a, &particle_b);
  if (distance_squared < 0.0)
    return -1;
  /* XXX: does this test NaN?? */
  if (distance_squared > 1000.0)
    return -2;

  /* check collision_time for parallel momenta => impossible collision */
  particle_a.set_momentum(0.1, 0.3, -0.1, 0.2);
  particle_b.set_momentum(0.1, 0.3, -0.1, 0.2);
  double time = collision_time(particle_a, particle_b);
  if (time >= 0.0)
    return -3;

  /* now check the Particles class itself */
  Particles particles;
  ParticleType piplus("pi+", 0.13957, -1.0, 211, 1, 1, 0);
  particles.add_type(piplus, 211);
  if (particles.types().size() != 1)
    return -4;
  size_t type_size = 0;
  for (std::map<int, ParticleType>::const_iterator
       i = particles.types().begin(); i != particles.types().end(); ++i) {
    printd("pdg %d mass: %g [GeV]\n", i->first, i->second.mass());
    type_size++;
  }
  if (type_size != 1)
    return -5;
  type_size = 0;
  ParticleType piminus("pi-", 0.13957, -1.0, -211, 1, -1, 0);
  for (std::map<int, ParticleType>::const_iterator
       i = particles.types().begin(); i != particles.types().end(); ++i) {
    printd("pdg %d mass: %g [GeV]\n", i->first, i->second.mass());
    type_size++;
  }
  if (type_size != 2)
    return -6;
  particles.add_type(piminus, -211);
  if (particles.types().size() != 2)
    return -7;
  particles.add_data(particle_a);
  if (particles.size() != 1)
    return -8;
  particles.add_data(particle_b);
  if (particles.size() != 2)
    return -9;
  double distance_squared_2 = particle_distance(particles.data_pointer(0),
    particles.data_pointer(1));
  if (distance_squared_2 < 0.0)
    return -10;
  if (distance_squared_2 - distance_squared < really_small)
    return -11;

  return 0;
}
