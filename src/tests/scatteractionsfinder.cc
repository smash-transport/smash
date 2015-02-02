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
#include "../include/action.h"
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

static void simple_propagate(Particles& particles, float dt) {
  for (ParticleData& data : particles.data()) {
    FourVector distance = FourVector(dt, data.velocity() * dt);
    data.set_4position(data.position() + distance);
  }
}

TEST(collision_order) {
  // create particles. The type doesn't matter at all because we will set a
  // different mass anyway and decays are switched off
  auto particle_a = create_smashon_particle(0);
  auto particle_b = create_smashon_particle(1);
  auto particle_c = create_smashon_particle(2);

  // set positions
  particle_a.set_4position(FourVector(0., 1., 0.1, 0.));
  particle_b.set_4position(FourVector(0., 0., 1., 0.));
  particle_c.set_4position(FourVector(0.,-1., 2., 0.));

  // set momenta, velocity should be very close to speed of light
  particle_a.set_4momentum(0.01, 0., 1., 0.);
  particle_b.set_4momentum(0.01, 1., 0., 0.);
  particle_c.set_4momentum(0.01, 1., 0., 0.);

  // put particles into list
  Particles particles;
  int id_a = particles.add_data(particle_a);
  int id_b = particles.add_data(particle_b);
  int id_c = particles.add_data(particle_c);

  // prepare scatteractionsfinder
  float elastic_parameter = 0.4;
  int testparticles = 1;
  ScatterActionsFinder finder(elastic_parameter, testparticles);

  // lists
  ParticleList search_list{particles.data().begin(), particles.data().end()};
  std::vector<const ParticleList*> neighbors_list; // empty for now

  // delta t (in fermi)
  float dt;

  // test for different times
  dt = 0.9;
  auto actions_1 = finder.find_possible_actions(search_list, neighbors_list, dt);
  COMPARE(actions_1.size(), 0) << "there shouldn't be a collision yet";

  dt = 1.0;
  auto actions_2 = finder.find_possible_actions(search_list, neighbors_list, dt);
  COMPARE(actions_2.size(), 1) << "there should already be a collision";

  dt = 2.0;
  auto actions_3 = finder.find_possible_actions(search_list, neighbors_list, dt);
  COMPARE(actions_3.size(), 2) << "there should be two collisions now";

  // perform actions from actions_3
  VERIFY(actions_3[0]->is_valid(particles)) << "first interaction is valid";
  size_t num_interactions = 0;
  actions_3[0]->perform(&particles, num_interactions);
  VERIFY(!actions_3[1]->is_valid(particles)) << "second interaction is not valid";
  COMPARE(num_interactions, 1)
      << "only one interaction should have been performed";

  // the position doesn't change from "perform", only the momenta
  //FourVector mom_a = particles.data(id_a).momentum();
}
