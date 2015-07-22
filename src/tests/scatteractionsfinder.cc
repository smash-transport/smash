/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include "setup.h"

#include <cstdio>

#include "../include/action.h"
#include "../include/constants.h"
#include "../include/particledata.h"
#include "../include/pdgcode.h"
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

  // delta t (in fermi)
  float dt;

  // test for different times
  dt = 0.9;
  auto actions_1 = finder.find_actions_in_cell(search_list, dt);
  COMPARE(actions_1.size(), 0u) << "timestep 0.9, expect no collision";

  dt = 1.0;
  auto actions_2 = finder.find_actions_in_cell(search_list, dt);
  COMPARE(actions_2.size(), 1u) << "timestep 1.0, expect 1 collision";

  dt = 2.0;
  auto actions_3 = finder.find_actions_in_cell(search_list, dt);
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

TEST(scatter_particle_pair_only_once) {
  // create two particles colliding head-on
  Particles p;
  p.insert(Test::smashon(Test::Momentum{0.11, 0., .1, 0.},
                         Test::Position{0., 1., .9, 1.}));
  p.insert(Test::smashon(Test::Momentum{0.11, 0., -.1, 0.},
                         Test::Position{0., 1., 1.1, 1.}));

  // prepare scatteractionsfinder
  const float radius = 0.11;                                        // in fm
  const float elastic_parameter = radius * radius * M_PI / fm2_mb;  // in mb
  const int testparticles = 1;
  ScatterActionsFinder finder(elastic_parameter, testparticles);
  ParticleList search_list = p.copy_to_vector();
  float dt = 0.9;  // fm/c

  // look for scatters, we expect one
  auto actions = finder.find_actions_in_cell(search_list, dt);
  COMPARE(actions.size(), 1u);

  // ok, the exepected Action exists, so let's perform it
  Action &action = *actions.front();
  action.generate_final_state();
  size_t processes = 0;
  action.perform(&p, processes);

  // because afterwards, the particles in p may not scatter again (they just
  // did)
  search_list = p.copy_to_vector();
  if (action.get_type() == ProcessType::Elastic) {
    // elastic scatters must keep the particle ids unchanged. Thus, the ids in
    // search_list must be 0 and 1.
    // TODO(weil): this test belongs into the Action unit tests, since this
    // tests the correct behavior of the Action class(es).
    COMPARE(action.outgoing_particles()[0].id(), 0);
    COMPARE(action.outgoing_particles()[1].id(), 1);
    COMPARE(search_list[0].id(), 0);
    COMPARE(search_list[1].id(), 1);
  }
  for (auto i = 10; i; --i) {  // make "sure" it's not the random numbers that
                               // supress the problem
    actions = finder.find_actions_in_cell(search_list, dt);
    COMPARE(actions.size(), 0u);
  }
}

TEST(find_next_action) {
  // let two particles collide head-on

  // x-position of the first particle (irrelevant for the calculation)
  constexpr double x_pos_1 = 0.0;
  // energy of the particles (irrelevant for the calculation)
  constexpr double energy = 1.0;
  // distance between the two particles
  constexpr double delta_x = 0.2;
  // velocity (in c) of the two particles
  constexpr double v = 0.5;

  // create particles
  Particles particles;
  particles.insert(Test::smashon(Test::Momentum{energy, energy * v, 0., 0.},
                                 Test::Position{0., x_pos_1, 1., 1.}));
  particles.insert(Test::smashon(Test::Momentum{energy, -energy*v, 0., 0.},
                                 Test::Position{0., x_pos_1+delta_x, 1., 1.}));


  // prepare scatteractionsfinder
  constexpr float radius = 0.11;                                        // in fm
  constexpr float elastic_parameter = radius * radius * M_PI / fm2_mb;  // in mb
  constexpr int testparticles = 1;
  ScatterActionsFinder finder(elastic_parameter, testparticles);

  // prepare list of particles that will be checked for possible actions
  ParticleList particle_list = particles.copy_to_vector();
  ActionList action_list = finder.find_actions_with_surrounding_particles(
      particle_list, particles, 10000.f);
  // we expect to find no actions because there are no surrounding particles
  COMPARE(action_list.size(), 0u);
  // remove one particle from the list so that the interaction can be found
  particle_list.pop_back();
  action_list = finder.find_actions_with_surrounding_particles(
      particle_list, particles, 10000.f);
  // we expect to find one collision between the two particles
  COMPARE(action_list.size(), 1u);
  ActionPtr action = std::move(action_list[0]);

  // calculate the expected time until the collision
  constexpr float collision_time = static_cast<float>(0.5*delta_x/v);
  // compare to what the action finder found
  FUZZY_COMPARE(action->time_of_execution(), collision_time);
}
