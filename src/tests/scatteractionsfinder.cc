/*
 *
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "setup.h"

#include <cstdio>

#include "../include/smash/action.h"
#include "../include/smash/constants.h"
#include "../include/smash/particledata.h"
#include "../include/smash/pdgcode.h"
#include "../include/smash/scatteractionsfinder.h"

using namespace smash;

TEST(init_particle_types) { Test::create_smashon_particletypes(); }

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
  particle_c.set_4position(FourVector(0., -1., 2., 0.));
  // particles d and e will also miss by 0.1 fm
  particle_d.set_4position(FourVector(0., 1.5, 0., 1.));
  particle_e.set_4position(FourVector(0., -0.1, 1.5, 1.));

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
  const double radius = 0.11;                                        // in fm
  const double elastic_parameter = radius * radius * M_PI / fm2_mb;  // in mb
  const std::vector<bool> has_interacted = {};
  ExperimentParameters exp_par = Test::default_parameters();
  Configuration config =
      Test::configuration("Collision_Term: {Elastic_Cross_Section: " +
                          std::to_string(elastic_parameter) + "}");
  ScatterActionsFinder finder(config, exp_par, has_interacted, 0, 0);

  // prepare lists
  ParticleList search_list = particles.copy_to_vector();

  // delta t (in fermi)
  double dt;

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
  uint32_t id_process = 1;

  // first action
  // verify that first action involves particle b
  const auto action_prtcls_1 = actions_3[0]->incoming_particles();
  VERIFY(std::find(action_prtcls_1.begin(), action_prtcls_1.end(),
                   particle_b) != action_prtcls_1.end());
  // check if action is valid
  VERIFY(actions_3[0]->is_valid(particles))
      << "expected: first interaction is valid";
  // perform action
  actions_3[0]->generate_final_state();
  actions_3[0]->perform(&particles, id_process);
  id_process++;

  // second action
  // verify that second action is *not* valid anymore
  VERIFY(!actions_3[1]->is_valid(particles))
      << "expected: second interaction is not valid";

  // third action
  // verify that third action involves particle d
  const auto action_prtcls_2 = actions_3[2]->incoming_particles();
  VERIFY(std::find(action_prtcls_2.begin(), action_prtcls_2.end(),
                   particle_d) != action_prtcls_2.end());
  // check if action is valid
  VERIFY(actions_3[2]->is_valid(particles))
      << "expected: third interaction is valid";
  // perform action
  actions_3[2]->generate_final_state();
  actions_3[2]->perform(&particles, id_process);
  id_process++;

  // final check
  COMPARE(id_process, 3u);
}

TEST(scatter_particle_pair_only_once) {
  // create two particles colliding head-on
  Particles p;
  p.insert(Test::smashon(Test::Momentum{0.11, 0., .1, 0.},
                         Test::Position{0., 1., .9, 1.}));
  p.insert(Test::smashon(Test::Momentum{0.11, 0., -.1, 0.},
                         Test::Position{0., 1., 1.1, 1.}));

  // prepare scatteractionsfinder
  const double radius = 0.11;                                        // in fm
  const double elastic_parameter = radius * radius * M_PI / fm2_mb;  // in mb
  const std::vector<bool> has_interacted = {};
  ExperimentParameters exp_par = Test::default_parameters();
  Configuration config =
      Test::configuration("Collision_Term: {Elastic_Cross_Section: " +
                          std::to_string(elastic_parameter) + "}");
  ScatterActionsFinder finder(config, exp_par, has_interacted, 0, 0);
  ParticleList search_list = p.copy_to_vector();
  double dt = 0.9;  // fm/c

  // look for scatters, we expect one
  auto actions = finder.find_actions_in_cell(search_list, dt);
  COMPARE(actions.size(), 1u);

  // ok, the exepected Action exists, so let's perform it
  Action &action = *actions.front();
  action.generate_final_state();
  const uint32_t id_process = 1;
  action.perform(&p, id_process);

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
                               // suppress the problem
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
  particles.insert(
      Test::smashon(Test::Momentum{energy, -energy * v, 0., 0.},
                    Test::Position{0., x_pos_1 + delta_x, 1., 1.}));

  // prepare scatteractionsfinder
  constexpr double radius = 0.11;  // in fm
  constexpr double elastic_parameter =
      radius * radius * M_PI / fm2_mb;  // in mb
  const std::vector<bool> has_interacted = {};
  ExperimentParameters exp_par = Test::default_parameters();
  Configuration config =
      Test::configuration("Collision_Term: {Elastic_Cross_Section: " +
                          std::to_string(elastic_parameter) + "}");
  ScatterActionsFinder finder(config, exp_par, has_interacted, 0, 0);

  // prepare list of particles that will be checked for possible actions
  ParticleList particle_list = particles.copy_to_vector();
  ActionList action_list = finder.find_actions_with_surrounding_particles(
      particle_list, particles, 10000.);
  // we expect to find no actions because there are no surrounding particles
  COMPARE(action_list.size(), 0u);
  // remove one particle from the list so that the interaction can be found
  particle_list.pop_back();
  action_list = finder.find_actions_with_surrounding_particles(
      particle_list, particles, 10000.);
  // we expect to find one collision between the two particles
  COMPARE(action_list.size(), 1u);
  ActionPtr action = std::move(action_list[0]);

  // calculate the expected time until the collision
  constexpr double collision_time = 0.5 * delta_x / v;
  // compare to what the action finder found
  FUZZY_COMPARE(action->time_of_execution(), collision_time);
}

TEST(increasing_scaling_factors) {
  constexpr double energy = 1.;
  constexpr double v = 0.5;
  constexpr double dx = 10;
  constexpr double dy = 1.;
  constexpr double time = -10.;
  // Cross section is 2 times area with a radius of impact the parameter dy
  constexpr double xsec = 4. * dy * dy * M_PI / fm2_mb;
  // Quadratic increase of scaling factor with time leads to linear increase of
  // the maximum distance, at which particles still collide, with time.
  constexpr double alpha = 2.;
  constexpr double delta_t_coll = dx / (2. * v);
  ParticleData p_a = ParticleData{ParticleType::find(0x661), 1};
  p_a.set_4position(FourVector(time, 0, 0, 0));
  p_a.set_4momentum(FourVector(energy, energy * v, 0, 0));
  ParticleData p_b = ParticleData{ParticleType::find(0x661), 2};
  p_b.set_4position(FourVector(time, dx, dy, 0));
  p_b.set_4momentum(FourVector(energy, -energy * v, 0, 0));
  // Set one particle to be half formed at time of collision.
  p_a.set_slow_formation_times(time, time + 2. * delta_t_coll);
  p_a.set_cross_section_scaling_factor(0.);
  // Set up scatter actions finder
  Configuration config = Test::configuration(
      "Collision_Term: {Elastic_Cross_Section: " + std::to_string(xsec) + "}");
  // For a power of larger than alpha, the particles should collide
  config.merge_yaml("Collision_Term: {Power_Particle_Formation: " +
                    std::to_string(alpha - 0.1) + "}");
  ExperimentParameters exp_par = Test::default_parameters();
  const std::vector<bool> has_interacted = {};
  ScatterActionsFinder finder(config, exp_par, has_interacted, 0, 0);
  COMPARE(finder.find_actions_in_cell({p_a, p_b}, 2. * delta_t_coll).size(),
          1u);
  // Set up a scatter actions finder so that the particles shouldn't collide
  Configuration config2 = Test::configuration(
      "Collision_Term: {Elastic_Cross_Section: " + std::to_string(xsec) + "}");
  config2.merge_yaml("Collision_Term: {Power_Particle_Formation: " +
                     std::to_string(alpha + 0.1) + "}");
  ScatterActionsFinder finder2(config2, exp_par, has_interacted, 0, 0);
  COMPARE(finder2.find_actions_in_cell({p_a, p_b}, 2. * delta_t_coll).size(),
          0u);
}
