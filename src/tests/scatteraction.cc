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
