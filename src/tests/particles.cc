/*
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */

#include "tests/unittest.h"

#include <cstdio>

#include "include/particles.h"
#include "include/constants.h"
#include "include/particledata.h"
#include "include/outputroutines.h"

TEST(everything) {
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
  VERIFY(!(distance_squared < 0.0));
  /* XXX: does this test NaN?? */
  VERIFY(!(distance_squared > 1000.0));

  /* check collision_time for parallel momenta => impossible collision */
  particle_a.set_momentum(0.1, 0.3, -0.1, 0.2);
  particle_b.set_momentum(0.1, 0.3, -0.1, 0.2);
  double time = collision_time(particle_a, particle_b);
  VERIFY(!(time >= 0.0));

  /* reset momenta for possible collision and compare to Particle class */
  particle_a.set_momentum(0.1, 10.0, 9.0, 8.0);
  particle_b.set_momentum(0.1, -10.0, -90.0, -80.0);
  distance_squared = particle_distance(&particle_a, &particle_b);
  VERIFY(!(distance_squared < 0.0));

  /* now check the Particles class itself */
  Particles particles;
  /* check addition of particle types */
  ParticleType piplus("pi+", 0.13957, -1.0, 211, 1, 1, 0);
  particles.add_type(piplus, 211);
  VERIFY(!(particles.types_size() != 1));
  size_t type_size = 0;
  for (std::map<int, ParticleType>::const_iterator
       i = particles.types_cbegin(); i != particles.types_cend(); ++i) {
    printd("pdg %d mass: %g [GeV]\n", i->first, i->second.mass());
    type_size++;
  }
  VERIFY(!(type_size != 1));
  type_size = 0;
  ParticleType piminus("pi-", 0.13957, -1.0, -211, 1, -1, 0);
  particles.add_type(piminus, -211);
  for (std::map<int, ParticleType>::const_iterator
       i = particles.types_cbegin(); i != particles.types_cend(); ++i) {
    printd("pdg %d mass: %g [GeV]\n", i->first, i->second.mass());
    type_size++;
  }
  VERIFY(!(type_size != 2));
  particles.add_type(piminus, -211);
  VERIFY(!(particles.types_size() != 2));

  /* check addition of particles */
  particles.add_data(particle_a);
  VERIFY(!(particles.size() != 1));
  particles.add_data(particle_b);
  VERIFY(!(particles.size() != 2));
  type_size = 0;
  for (std::map<int, ParticleData>::const_iterator
       i = particles.cbegin(); i != particles.cend(); ++i) {
    printd("id %d: pdg %d\n", i->first, i->second.pdgcode());
    /* check that id and and own id are the same */
    VERIFY(!(i->first != i->second.id()));
    type_size++;
  }
  VERIFY(!(type_size != 2));
  VERIFY(!(particles.empty()));
  VERIFY(!(particles.count(1) != 1));

  /* check usage particle data */
  double distance_squared_2 = particle_distance(particles.data_pointer(0),
    particles.data_pointer(1));
  printd("%g versus %g\n", distance_squared, distance_squared_2);
  VERIFY(!(distance_squared_2 < 0.0));
  VERIFY(!(distance_squared_2 - distance_squared > really_small));
}
