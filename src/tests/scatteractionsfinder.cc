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

TEST(collision_order) {
  // create particles, the type doesn't matter at all because we will set a
  // different mass anyway and decays are switched off
  auto particle_a = create_smashon_particle(0);
  auto particle_b = create_smashon_particle(1);
  auto particle_c = create_smashon_particle(2);
  auto particle_d = create_smashon_particle(3);
  auto particle_e = create_smashon_particle(4);

  // set positions
  // particle a is set such that it will miss particles b and c by 0.1 fm
  particle_a.set_4position(FourVector(0., 1., 0.1, 0.));
  particle_b.set_4position(FourVector(0., 0., 1., 0.));
  particle_c.set_4position(FourVector(0.,-1., 2., 0.));
  // particles d and e will also miss by 0.1 fm
  particle_d.set_4position(FourVector(0., 1.5, 0., 1.));
  particle_e.set_4position(FourVector(0.,-0.1, 1.5, 1.));

  // set momenta
  // velocity should be very close to speed of light to make calculations easier
  particle_a.set_4momentum(0.01, 0., 1., 0.);
  particle_b.set_4momentum(0.01, 1., 0., 0.);
  particle_c.set_4momentum(0.01, 1., 0., 0.);
  particle_d.set_4momentum(0.01, 0., 1., 0.);
  particle_e.set_4momentum(0.01, 1., 0., 0.);

  // put particles into list
  Particles particles;
  particles.insert(particle_a);
  particles.insert(particle_b);
  particles.insert(particle_c);
  particles.insert(particle_d);
  particles.insert(particle_e);

  // prepare scatteractionsfinder
  const float radius = 0.11; // in fm
  const float elastic_parameter = radius*radius*M_PI/fm2_mb; // in mb
  const int testparticles = 1;
  ScatterActionsFinder finder(elastic_parameter, testparticles);

  // prepare lists
  ParticleList search_list = particles.copy_to_vector();
  std::vector<const ParticleList*> neighbors_list; // empty for now

  // delta t (in fermi)
  float dt;

  // test for different times
  dt = 0.9;
  auto actions_1 = finder.find_possible_actions(search_list,neighbors_list, dt);
  COMPARE(actions_1.size(), 0u) << "timestep 0.9, expect no collision";

  dt = 1.0;
  auto actions_2 = finder.find_possible_actions(search_list,neighbors_list, dt);
  COMPARE(actions_2.size(), 1u) << "timestep 1.0, expect 1 collision";

  dt = 2.0;
  auto actions_3 = finder.find_possible_actions(search_list,neighbors_list, dt);
  COMPARE(actions_3.size(), 3u) << "timestep 2.0, expect 3 collisions";

  // perform actions from actions_3

  size_t num_interactions = 0;

  // first action
  // verify that first action involves particle b
  const auto action_prtcls_1 = actions_3[0]->incoming_particles();
  VERIFY(std::find(action_prtcls_1.begin(), action_prtcls_1.end(), particle_b)
          != action_prtcls_1.end());
  // check if action is valid
  VERIFY(actions_3[0]->is_valid(particles))
      << "expected: first interaction is valid";
  // perform action
  actions_3[0]->generate_final_state();
  actions_3[0]->perform(&particles, num_interactions);

  // second action
  // verify that second action is *not* valid anymore
  VERIFY(!actions_3[1]->is_valid(particles))
      << "expected: second interaction is not valid";

  // third action
  // verify that third action involves particle d
  const auto action_prtcls_2 = actions_3[2]->incoming_particles();
  VERIFY(std::find(action_prtcls_2.begin(), action_prtcls_2.end(), particle_d)
          != action_prtcls_2.end());
  // check if action is valid
  VERIFY(actions_3[2]->is_valid(particles))
      << "expected: third interaction is valid";
  // perform action
  actions_3[2]->generate_final_state();
  actions_3[2]->perform(&particles, num_interactions);

  // final check
  COMPARE(num_interactions, 2u);
}
