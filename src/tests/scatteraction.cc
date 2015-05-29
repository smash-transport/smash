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

#include "../include/scatteraction.h"
#include "../include/scatteractionbaryonmeson.h"

using namespace Smash;
using Smash::Test::Position;
using Smash::Test::Momentum;

TEST(init_particle_types) {
  Test::create_smashon_particletypes();
}

TEST(sorting) {
  const auto a = Test::smashon(Position{0., -0.1, 0., 0.},
                               Momentum{1.1, 1.0, 0., 0.});
  const auto b = Test::smashon(Position{0., 0.1, 0., 0.},
                               Momentum{1.1, 1.0, 0., 0.});
  constexpr float time1 = 1.f;
  ScatterAction act1(a, b, time1);
  COMPARE(act1.get_interaction_point(), FourVector(0., 0., 0., 0.));

  constexpr float time2 = 1.1f;
  ScatterAction act2(a, b, time2);
  VERIFY(act1 < act2);
}

TEST(elastic_collision) {
  // put particles in list
  Particles particles;
  particles.insert(Test::smashon(Position{0., -0.1, 0., 0.},
                               Momentum{1.1, 1.0, 0., 0.}));
  particles.insert(Test::smashon(Position{0., 0.1, 0., 0.},
                               Momentum{1.1, 1.0, 0., 0.}));
  const ParticleData& a = particles.front();
  const ParticleData& b = particles.back();

  // create action
  constexpr float time = 1.f;
  ScatterAction act(a, b, time);
  VERIFY(act.is_valid(particles));

  // add elastic channel
  constexpr float sigma = 10.0;
  act.add_collision(act.elastic_cross_section(sigma));

  // check cross section
  COMPARE(act.cross_section(), sigma);

  // generate final state
  act.generate_final_state();

  // verify that the action is indeed elastic
  COMPARE(act.get_type(), ProcessType::Elastic);

  // verify that particles didn't change in the collision
  ParticleList in = act.incoming_particles();
  const ParticleList& out = act.outgoing_particles();
  VERIFY((in[0] == out[0] && in[1] == out[1])
         || (in[0] == out[1] && in[1] == out[0]));

  // perform the action
  size_t id_process = 1u;
  act.perform(&particles, id_process);

  // action should not be valid anymore
  VERIFY(!act.is_valid(particles));

  // verify that the particles don't change in the particle list
  VERIFY((in[0] == particles.front() && in[1] == particles.back())
         || (in[0] == particles.back() && in[1] == particles.front()));
}

TEST(outgoing_valid) {
  Test::create_actual_particletypes();
  Test::create_actual_decaymodes();
  // create a proton and a pion
  ParticleData p1{ParticleType::find(0x2212)};
  ParticleData p2{ParticleType::find(0x111)};
  // set position
  p1.set_4position(Position{0., -0.1, 0., 0.});
  p2.set_4position(Position{0., 0.1, 0., 0.});
  // set momenta
  constexpr double p_x = 0.1;
  const double mass = p1.pole_mass();
  const double energy = std::sqrt(mass*mass+p_x*p_x);
  p1.set_4momentum(Momentum{energy, p_x, 0., 0.});
  p2.set_4momentum(Momentum{energy, -p_x, 0., 0.});

  // put in particles object
  Particles particles;
  particles.insert(p1);
  particles.insert(p2);

  // get valid copies back
  ParticleList plist = particles.copy_to_vector();
  auto p1_copy = plist[0];
  auto p2_copy = plist[1];
  VERIFY(particles.is_valid(p1_copy) && particles.is_valid(p2_copy));

  // construct action
  ScatterActionPtr act;
  act = make_unique<ScatterActionBaryonMeson>(p1_copy, p2_copy, 0.2f);
  VERIFY(act != nullptr);
  COMPARE(p2_copy.type(), ParticleType::find(0x111));

  // add processes
  constexpr float elastic_parameter = 0.f;  // don't include elastic scattering
  act->add_all_processes(elastic_parameter);
  VERIFY(act->cross_section() > 0.f);

  // perform actions
  VERIFY(act->is_valid(particles));
  act->generate_final_state();
  VERIFY(act->get_type() != ProcessType::Elastic);
  size_t id_process = 0u;
  act->perform(&particles, id_process);
  COMPARE(id_process, 1u);

  // check the outgoing particles
  const ParticleList& outgoing_particles = act->outgoing_particles();
  VERIFY(outgoing_particles.size() > 0u);  // should be at least one
  VERIFY(particles.is_valid(outgoing_particles[0]));
}
