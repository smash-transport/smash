/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"

#include <cstdio>

#include "../include/particledata.h"
#include "../include/pdgcode.h"
#include "../include/scatteraction.h"
#include "../include/scatteractionsfinder.h"

using namespace Smash;

TEST(init_particle_types) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PDG\n"
      "smashon 0.123 1.2 661\n");
}

static ParticleData create_smashon_particle(int id = -1) {
  return ParticleData{ParticleType::find(0x661), id};
}

TEST(everything) {
  ParticleData particle_a = create_smashon_particle(0),
               particle_b = create_smashon_particle(1);
  particle_a.set_4position(FourVector(1., 1., 1., 1.));
  particle_b.set_4position(FourVector(2., 2., 2., 2.));

  /* test collision_time:
   * parallel momenta => impossible collision */
  particle_a.set_4momentum(0.1, 0.3, -0.1, 0.2);
  particle_b.set_4momentum(0.1, 0.3, -0.1, 0.2);
  double time = ScatterActionsFinder::collision_time(particle_a, particle_b);
  VERIFY(time < 0.0);

  /* test particle_distance:
   * particles with null momenta */
  particle_a.set_4momentum(0.1, 0.0, 0.0, 0.0);
  particle_b.set_4momentum(0.1, 0.0, 0.0, 0.0);
  ScatterAction *act = new ScatterAction(particle_a, particle_b, 0.);
  double distance_squared = act->particle_distance();
  VERIFY(distance_squared >= 0.);
  VERIFY(distance_squared <= 100.);
  delete(act);

  /* particles with finite momenta */
  particle_a.set_4momentum(0.1, 10.0, 9.0, 8.0);
  particle_b.set_4momentum(0.1, -10.0, -90.0, -80.0);
  act = new ScatterAction(particle_a, particle_b, 0.);
  distance_squared = act->particle_distance();
  VERIFY(distance_squared >= 0.);
  delete act;
}
