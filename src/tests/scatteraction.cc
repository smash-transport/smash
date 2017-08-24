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

#include "../include/scatteractionbaryonmeson.h"
#include "../include/scatteractionbaryonbaryon.h"

using namespace Smash;
using Smash::Test::Position;
using Smash::Test::Momentum;

TEST(init_particle_types) {
  Test::create_actual_particletypes();
  Test::create_actual_decaymodes();
}

constexpr double r_x = 0.1;
const FourVector pos_a = Position{0., -r_x, 0., 0.};
const FourVector pos_b = Position{0.,  r_x, 0., 0.};
const FourVector middle = (pos_a + pos_b) / 2.;

TEST(sorting) {
  ParticleData a{ParticleType::find(0x111)};  // pi0
  a.set_4position(pos_a);
  a.set_4momentum(Momentum{1.1, 1.0, 0., 0.});

  ParticleData b{ParticleType::find(0x111)};  // pi0
  b.set_4position(pos_b);
  a.set_4momentum(Momentum{1.1, 1.0, 0., 0.});

  constexpr double time1 = 1.;
  ScatterAction act1(a, b, time1);
  COMPARE(act1.get_interaction_point(), middle);

  constexpr double time2 = 1.1;
  ScatterAction act2(a, b, time2);
  VERIFY(act1 < act2);
}

TEST(elastic_collision) {
  // put particles in list
  Particles particles;
  ParticleData a{ParticleType::find(0x211)};  // pi+
  a.set_4position(pos_a);
  a.set_4momentum(Momentum{1.1, 1.0, 0., 0.});
  a.set_history(3, 1, ProcessType::None, 1.2, ParticleList{});

  ParticleData b{ParticleType::find(0x211)};  // pi+
  b.set_4position(pos_b);
  b.set_4momentum(Momentum{1.1, 1.0, 0., 0.});
  b.set_history(3, 1, ProcessType::None, 1.2, ParticleList{});

  a = particles.insert(a);
  b = particles.insert(b);

  // create action
  constexpr double time = 1.;
  ScatterAction act(a, b, time);
  ScatterAction act_copy(a, b, time);
  VERIFY(act.is_valid(particles));
  VERIFY(act_copy.is_valid(particles));

  // add elastic channel
  constexpr double sigma = 10.0;
  constexpr bool strings_switch = false;
  constexpr NNbarTreatment nnbar_treatment = NNbarTreatment::NoAnnihilation;
  act.add_all_processes(sigma, true, 0., strings_switch, nnbar_treatment);

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

  // verify that the particles keep their positions after elastic scattering
  COMPARE(out[0].position(), pos_a);
  COMPARE(out[1].position(), pos_b);

  // perform the action
  COMPARE(particles.front().id_process(), 1u);
  uint32_t id_process = 2;
  act.perform(&particles, id_process);
  id_process++;
  // check id_process
  COMPARE(particles.front().id_process(), 2u);
  COMPARE(particles.back().id_process(), 2u);
  COMPARE(id_process, 3u);

  // action should not be valid anymore
  VERIFY(!act.is_valid(particles));
  VERIFY(!act_copy.is_valid(particles));

  // verify that the particles don't change in the particle list
  VERIFY((in[0] == particles.front() && in[1] == particles.back())
         || (in[0] == particles.back() && in[1] == particles.front()));
}

TEST(outgoing_valid) {
  // create a proton and a pion
  ParticleData p1{ParticleType::find(0x2212)};
  ParticleData p2{ParticleType::find(0x111)};
  // set position
  p1.set_4position(pos_a);
  p2.set_4position(pos_b);
  // set momenta
  constexpr double p_x = 0.1;
  p1.set_4momentum(p1.pole_mass(), p_x, 0., 0.);
  p2.set_4momentum(p2.pole_mass(), -p_x, 0., 0.);

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
  act = make_unique<ScatterActionBaryonMeson>(p1_copy, p2_copy, 0.2);
  VERIFY(act != nullptr);
  COMPARE(p2_copy.type(), ParticleType::find(0x111));

  // add processes
  constexpr double elastic_parameter = 0.;  // don't include elastic scattering
  constexpr bool strings_switch = false;
  constexpr NNbarTreatment nnbar_treatment = NNbarTreatment::NoAnnihilation;
  act->add_all_processes(elastic_parameter, true, 0., strings_switch,
                         nnbar_treatment);

  VERIFY(act->cross_section() > 0.);

  // perform actions
  VERIFY(act->is_valid(particles));
  act->generate_final_state();
  VERIFY(act->get_type() != ProcessType::Elastic);
  const uint32_t id_process = 1;
  act->perform(&particles, id_process);
  COMPARE(id_process, 1u);

  // check the outgoing particles
  const ParticleList& outgoing_particles = act->outgoing_particles();
  VERIFY(outgoing_particles.size() > 0u);  // should be at least one
  VERIFY(particles.is_valid(outgoing_particles[0]));
  VERIFY(outgoing_particles[0].id() > p1_copy.id());
  VERIFY(outgoing_particles[0].id() > p2_copy.id());
  // verify that particle is placed in the middle between the incoming ones
  COMPARE(outgoing_particles[0].position(), middle);
}

TEST(pythia_running) {
  // create two protons
  ParticleData p1{ParticleType::find(0x2212)};
  ParticleData p2{ParticleType::find(0x2212)};
  // set position
  p1.set_4position(pos_a);
  p2.set_4position(pos_b);
  // set momenta
  constexpr double p_x = 3.0;
  p1.set_4momentum(p1.pole_mass(), p_x, 0., 0.);
  p2.set_4momentum(p2.pole_mass(), -p_x, 0., 0.);

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
  std::set<IncludedReactions> incl_2to2 = {};
  act = make_unique<ScatterActionBaryonBaryon>(p1_copy, p2_copy, 0.2, false,
                                               incl_2to2, 1.0);
  VERIFY(act != nullptr);
  COMPARE(p2_copy.type(), ParticleType::find(0x2212));

  // add processes
  constexpr double elastic_parameter = 0.;  // don't include elastic scattering
  constexpr bool strings_switch = true;
  constexpr NNbarTreatment nnbar_treatment = NNbarTreatment::NoAnnihilation;
  act->add_all_processes(elastic_parameter, false, 0., strings_switch,
                         nnbar_treatment);

  VERIFY(act->cross_section() > 0.);

  // perform actions
  VERIFY(act->is_valid(particles));
  act->generate_final_state();
  VERIFY(act->get_type() != ProcessType::Elastic);
  VERIFY(act->get_type() == ProcessType::String);
  const uint32_t id_process = 1;
  act->perform(&particles, id_process);
  COMPARE(id_process, 1u);

  // check the outgoing particles
  const ParticleList& outgoing_particles = act->outgoing_particles();
  VERIFY(outgoing_particles.size() > 0u);  // should be at least one
  VERIFY(particles.is_valid(outgoing_particles[0]));
  VERIFY(outgoing_particles[0].id() > p1_copy.id());
  VERIFY(outgoing_particles[0].id() > p2_copy.id());
}

TEST(update_incoming) {
  // put particles in list
  Particles particles;
  ParticleData a{ParticleType::find(0x211)};  // pi+
  a.set_4position(pos_a);
  a.set_4momentum(Momentum{1.1, 1.0, 0., 0.});

  ParticleData b{ParticleType::find(0x211)};  // pi+
  b.set_4position(pos_b);
  b.set_4momentum(Momentum{1.1, -1.0, 0., 0.});

  a = particles.insert(a);
  b = particles.insert(b);

  // create action
  constexpr double time = 0.2;
  ScatterAction act(a, b, time);
  VERIFY(act.is_valid(particles));

  // add elastic channel
  constexpr double sigma = 10.0;
  bool string_switch = true;
  NNbarTreatment nnbar_treatment = NNbarTreatment::NoAnnihilation;
  act.add_all_processes(sigma, true, 0., string_switch, nnbar_treatment);

  // change the position of one of the particles
  const FourVector new_position(0.1, 0., 0., 0.);
  particles.front().set_4position(new_position);

  // update the action
  act.update_incoming(particles);
  COMPARE(act.incoming_particles()[0].position(), new_position);
}
