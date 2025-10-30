/*
 *
 *    Copyright (c) 2015-2020,2022,2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/scatteractionsfinder.h"

#include <cstdio>

#include "setup.h"
#include "smash/action.h"
#include "smash/constants.h"
#include "smash/particledata.h"
#include "smash/pdgcode.h"

using namespace smash;

TEST(init_particle_types) { Test::create_smashon_particletypes(); }

static ParticleData create_smashon_particle(int id = -1) {
  return ParticleData{ParticleType::find(0x661), id};
}

static Configuration create_configuration_for_tests(double cross_section) {
  const std::string tmp{R"(
    Collision_Term:
      Elastic_Cross_Section: )" +
                        std::to_string(cross_section)};
  return Configuration{tmp.c_str()};
}

TEST(collision_order) {
  // create particles, the type doesn't matter at all because we will set a
  // different mass anyway and decays are switched off
  auto particle_a = create_smashon_particle(0);
  auto particle_b = create_smashon_particle(1);
  auto particle_c = create_smashon_particle(2);
  auto particle_d = create_smashon_particle(3);
  auto particle_e = create_smashon_particle(4);

  // set formation times
  particle_a.set_formation_time(0.);
  particle_b.set_formation_time(0.);
  particle_c.set_formation_time(0.);
  particle_d.set_formation_time(0.);
  particle_e.set_formation_time(0.);

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
  ExperimentParameters exp_par = Test::default_parameters();
  Configuration config = create_configuration_for_tests(elastic_parameter);
  ScatterActionsFinder finder(config, exp_par);

  // prepare lists
  ParticleList search_list = particles.copy_to_vector();

  // no grid means no grid cell volume
  const double grid_cell_vol = 0.0;

  // frozen Fermi motion is not tested, so we do not need the beam momentum
  const std::vector<FourVector> beam_mom = {};

  // delta t (in fermi)
  double dt;

  // test for different times
  dt = 0.9;
  auto actions_1 =
      finder.find_actions_in_cell(search_list, dt, grid_cell_vol, beam_mom);
  COMPARE(actions_1.size(), 0u) << "timestep 0.9, expect no collision";

  dt = 1.0;
  auto actions_2 =
      finder.find_actions_in_cell(search_list, dt, grid_cell_vol, beam_mom);
  COMPARE(actions_2.size(), 1u) << "timestep 1.0, expect 1 collision";

  dt = 2.0;
  auto actions_3 =
      finder.find_actions_in_cell(search_list, dt, grid_cell_vol, beam_mom);
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
  ExperimentParameters exp_par = Test::default_parameters();
  Configuration config = create_configuration_for_tests(elastic_parameter);
  ScatterActionsFinder finder(config, exp_par);
  ParticleList search_list = p.copy_to_vector();
  double dt = 0.9;                   // fm
  const double grid_cell_vol = 0.0;  // no grid

  // look for scatters, we expect one
  auto actions =
      finder.find_actions_in_cell(search_list, dt, grid_cell_vol, {});
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
    actions = finder.find_actions_in_cell(search_list, dt, grid_cell_vol, {});
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
  ExperimentParameters exp_par = Test::default_parameters();
  Configuration config = create_configuration_for_tests(elastic_parameter);
  ScatterActionsFinder finder(config, exp_par);

  // prepare list of particles that will be checked for possible actions
  ParticleList particle_list = particles.copy_to_vector();
  ActionList action_list = finder.find_actions_with_surrounding_particles(
      particle_list, particles, 10000., {});
  // we expect to find no actions because there are no surrounding particles
  COMPARE(action_list.size(), 0u);
  // remove one particle from the list so that the interaction can be found
  particle_list.pop_back();
  action_list = finder.find_actions_with_surrounding_particles(
      particle_list, particles, 10000., {});
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
  constexpr double grid_cell_vol = 0.0;  // no grid
  constexpr double delta_t_coll = dx / (2. * v);
  ParticleData p_a = ParticleData{ParticleType::find(0x661), 1};
  p_a.set_4position(FourVector(time, 0, 0, 0));
  p_a.set_4momentum(FourVector(energy, energy * v, 0, 0));
  ParticleData p_b = ParticleData{ParticleType::find(0x661), 2};
  p_b.set_4position(FourVector(time, dx, dy, 0));
  p_b.set_4momentum(FourVector(energy, -energy * v, 0, 0));
  // Set one particle to be half formed at time of collision.
  p_a.set_slow_formation_times(time, time + 2. * delta_t_coll);
  p_b.set_formation_time(time);
  p_a.set_cross_section_scaling_factor(0.);
  // Set up scatter actions finder
  Configuration config = create_configuration_for_tests(xsec);
  // For a power of larger than alpha, the particles should collide.
  ParticleData::formation_power_ = alpha - 0.1;
  ExperimentParameters exp_par = Test::default_parameters();
  ScatterActionsFinder finder(config, exp_par);
  COMPARE(finder
              .find_actions_in_cell({p_a, p_b}, 2. * delta_t_coll,
                                    grid_cell_vol, {})
              .size(),
          1u);
  // For a Power smaller than alpha, the particles should not collide.
  ParticleData::formation_power_ = alpha + 0.1;
  COMPARE(finder
              .find_actions_in_cell({p_a, p_b}, 2. * delta_t_coll,
                                    grid_cell_vol, {})
              .size(),
          0u);
}

TEST(check_stochastic_collision) {
  // create two particles with some momentum
  Particles p;
  p.insert(Test::smashon(Test::Momentum{0.11, 0., .1, 0.},
                         Test::Position{0., 1., .9, 1.}));
  p.insert(Test::smashon(Test::Momentum{0.11, 0., -.1, 0.},
                         Test::Position{0., 1., 1.1, 1.}));

  // use fictious numbers for collisions search
  const double grid_cell_vol = 8.0;
  const double dt = 0.1;
  const int testparticles = 1;

  // prepare scatteractionsfinder
  const double elastic_parameter = 10.0;  // in mb
  ExperimentParameters exp_par = Test::default_parameters(
      testparticles, dt, CollisionCriterion::Stochastic);
  Configuration config = create_configuration_for_tests(elastic_parameter);
  ScatterActionsFinder finder(config, exp_par);

  // prepare lists
  ParticleList search_list = p.copy_to_vector();

  // calcluate the relative velocity by hand
  const double m1 = search_list[0].effective_mass();
  const double m2 = search_list[1].effective_mass();
  FourVector mom(0.0, 0.0, 0.0, 0.0);
  mom += search_list[0].momentum();
  mom += search_list[1].momentum();
  const double m_s = mom.sqr();
  const double lamb = Action::lambda_tilde(m_s, m1 * m1, m2 * m2);
  const double v_rel = std::sqrt(lamb) / (2. * search_list[0].momentum().x0() *
                                          search_list[1].momentum().x0());

  const int N_samples = 1E6;
  int found_actions = 0;
  for (int i = 0; i < N_samples; i++) {
    auto actions =
        finder.find_actions_in_cell(search_list, dt, grid_cell_vol, {});
    found_actions += actions.size();
  }

  // probability of finding an action
  const double ratio_found =
      static_cast<double>(found_actions) / static_cast<double>(N_samples);

  // calculate probability of stochastic criterion by hand
  const double prob = elastic_parameter * fm2_mb * v_rel * dt / grid_cell_vol;

  // compare probability to the probability of finding an action
  COMPARE_RELATIVE_ERROR(ratio_found, prob, 0.05);
}
