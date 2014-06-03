/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include "../include/modusdefault.h"
#include "../include/experiment.h"

using namespace Smash;

TEST(sanity) {
  ModusDefault m;
  Particles * P;
  COMPARE(m.sanity_check(P), 0);
}

TEST(propagate) {
  ModusDefault m;
  Particles P{"", ""};
  P.reset();
  // particle_static shouldn't move:
  ParticleData particle_static;
  FourVector mom_1(4.0, 0.0, 0.0, 0.0);
  FourVector pos_1(5.0, 6.0, 6.7, 8.9);
  particle_static.set_momentum(mom_1);
  particle_static.set_position(pos_1);
  // particle_light should move with speed of light:
  ParticleData particle_light;
  FourVector mom_2(sqrt(.02), 0.1, -.1, 0.0);
  COMPARE_ABSOLUTE_ERROR(mom_2.sqr(), 0.0, 1e-17);
  FourVector pos_2(5.0, 4.3, 2.1, -.1);
  particle_light.set_momentum(mom_2);
  particle_light.set_position(pos_2);
  ParticleData particle_slow;
  FourVector mom_3(sqrt(1.0 + .01 + .04 + .09), 0.1, 0.2, 0.3);
  FourVector pos_3(8.7, 6.5, 4.3, 2.1);
  particle_slow.set_momentum(mom_3);
  particle_slow.set_position(pos_3);
  int id_static = P.add_data(particle_static);
  int id_light  = P.add_data(particle_light);
  int id_slow   = P.add_data(particle_slow);
  OutputsList out;
  // clock, output interval, cross-section, testparticles
  ExperimentParameters param{{0.f, 1.f}, 1.f, 0.0, 1};
  m.propagate(&P, param, out);
  // after propagation: Momentum should be unchanged.
  COMPARE(P.data(id_static).momentum().x0(), mom_1.x0());
  COMPARE(P.data(id_static).momentum().x1(), mom_1.x1());
  COMPARE(P.data(id_static).momentum().x2(), mom_1.x2());
  COMPARE(P.data(id_static).momentum().x3(), mom_1.x3());
  COMPARE(P.data(id_light ).momentum().x0(), mom_2.x0());
  COMPARE(P.data(id_light ).momentum().x1(), mom_2.x1());
  COMPARE(P.data(id_light ).momentum().x2(), mom_2.x2());
  COMPARE(P.data(id_light ).momentum().x3(), mom_2.x3());
  COMPARE(P.data(id_slow  ).momentum().x0(), mom_3.x0());
  COMPARE(P.data(id_slow  ).momentum().x1(), mom_3.x1());
  COMPARE(P.data(id_slow  ).momentum().x2(), mom_3.x2());
  COMPARE(P.data(id_slow  ).momentum().x3(), mom_3.x3());
  // position should be updated:
  double time = param.new_particle_time();
  double distance_2 = param.timestep_duration() / std::sqrt(2.0);
  ThreeVector distance_3 = mom_3.threevec() / mom_3.x0();
  FUZZY_COMPARE(P.data(id_static).position().x0(), time);
  FUZZY_COMPARE(P.data(id_static).position().x1(), pos_1.x1());
  FUZZY_COMPARE(P.data(id_static).position().x2(), pos_1.x2());
  FUZZY_COMPARE(P.data(id_static).position().x3(), pos_1.x3());
  FUZZY_COMPARE(P.data(id_light ).position().x0(), time);
  FUZZY_COMPARE(P.data(id_light ).position().x1(), pos_2.x1() + distance_2);
  FUZZY_COMPARE(P.data(id_light ).position().x2(), pos_2.x2() - distance_2);
  FUZZY_COMPARE(P.data(id_light ).position().x3(), pos_2.x3());
  FUZZY_COMPARE(P.data(id_slow  ).position().x0(), time);
  FUZZY_COMPARE(P.data(id_slow  ).position().x1(), pos_3.x1() + distance_3.x1());
  FUZZY_COMPARE(P.data(id_slow  ).position().x2(), pos_3.x2() + distance_3.x2());
  FUZZY_COMPARE(P.data(id_slow  ).position().x3(), pos_3.x3() + distance_3.x3());
}
