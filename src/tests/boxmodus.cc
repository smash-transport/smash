/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include "../include/boxmodus.h"
#include "../include/experiment.h"
#include "../include/configuration.h"

#include <boost/filesystem.hpp>

using namespace Smash;

TEST(sanity) {
  Configuration conf(TEST_CONFIG_PATH);
  conf["Modi"]["Box"]["INITIAL_CONDITION"] = 1;
  conf["Modi"]["Box"]["LENGTH"] = 1.0;
  conf["Modi"]["Box"]["TEMPERATURE"] = 0.13;
  conf["Modi"]["Box"]["START_TIME"] = 0.2;
  ExperimentParameters param{{0.f, 1.f}, 1.f, 0.0, 1};
  BoxModus box(conf["Modi"], param);

  Particles P{"", ""};
  ParticleData particle_outside;
  FourVector pos_1(5.0, 1.2, -.1, 0.2);
  particle_outside.set_position(pos_1);
  int id1 = P.add_data(particle_outside);
  ParticleData particle_inside;
  FourVector pos_2(7.0, 0.8, 0.7, 0.9);
  particle_inside.set_position(pos_2);
  int id2 = P.add_data(particle_inside);

  UnitTest::setFuzzyness<double>(2);
  COMPARE(box.sanity_check(&P), 1);
  COMPARE(P.data(id1).position().x0(), pos_1.x0());
  COMPARE(P.data(id1).position().x1(), pos_1.x1() - 1.0);
  COMPARE(P.data(id1).position().x2(), pos_1.x2() + 1.0);
  COMPARE(P.data(id1).position().x3(), pos_1.x3());
  COMPARE(P.data(id2).position().x0(), pos_2.x0());
  COMPARE(P.data(id2).position().x1(), pos_2.x1());
  COMPARE(P.data(id2).position().x2(), pos_2.x2());
  COMPARE(P.data(id2).position().x3(), pos_2.x3());
}

TEST(propagate_as_default_does) {
  Configuration conf(TEST_CONFIG_PATH);
  conf["Modi"]["Box"]["INITIAL_CONDITION"] = 1;
  conf["Modi"]["Box"]["LENGTH"] = 1.0;
  conf["Modi"]["Box"]["TEMPERATURE"] = 0.13;
  conf["Modi"]["Box"]["START_TIME"] = 0.2;
  ExperimentParameters param{{0.f, 0.1f}, 1.f, 0.0, 1};
  BoxModus box(conf["Modi"], param);

  ModusDefault m;
  // one particle list per Modus.
  Particles Pbox{"", ""};
  Particles Pdef{"", ""};

  // we check first if propagate does the same as for ModusDefault where
  // it should:

  // particle_static shouldn't move:
  ParticleData particle_static;
  FourVector mom_1(4.0, 0.0, 0.0, 0.0);
  FourVector pos_1(0.0, 0.6, 0.7, 0.9);
  particle_static.set_momentum(mom_1);
  particle_static.set_position(pos_1);
  // particle_light should move with speed of light:
  ParticleData particle_light;
  FourVector mom_2(sqrt(.02), 0.1, -.1, 0.0);
  COMPARE_ABSOLUTE_ERROR(mom_2.sqr(), 0.0, 1e-17);
  FourVector pos_2(5.0, 0.3, 0.1, .5);
  particle_light.set_momentum(mom_2);
  particle_light.set_position(pos_2);
  ParticleData particle_slow;
  FourVector mom_3(sqrt(1.0 + .01 + .04 + .09), 0.1, 0.2, 0.3);
  FourVector pos_3(8.7, 0.1, 0.1, 0.1);
  particle_slow.set_momentum(mom_3);
  particle_slow.set_position(pos_3);
  int box_static = Pbox.add_data(particle_static);
  int box_light  = Pbox.add_data(particle_light);
  int box_slow   = Pbox.add_data(particle_slow);
  int def_static = Pdef.add_data(particle_static);
  int def_light  = Pdef.add_data(particle_light);
  int def_slow   = Pdef.add_data(particle_slow);
  OutputsList out;
    m.propagate(&Pdef, param, out);
  box.propagate(&Pbox, param, out);
  // for those particles, Box and Default Modi should propagate equally:
  COMPARE(Pbox.data(box_static).momentum().x0(), Pdef.data(def_static).momentum().x0());
  COMPARE(Pbox.data(box_static).momentum().x1(), Pdef.data(def_static).momentum().x1());
  COMPARE(Pbox.data(box_static).momentum().x2(), Pdef.data(def_static).momentum().x2());
  COMPARE(Pbox.data(box_static).momentum().x3(), Pdef.data(def_static).momentum().x3());
  COMPARE(Pbox.data(box_light ).momentum().x0(), Pdef.data(def_light ).momentum().x0());
  COMPARE(Pbox.data(box_light ).momentum().x1(), Pdef.data(def_light ).momentum().x1());
  COMPARE(Pbox.data(box_light ).momentum().x2(), Pdef.data(def_light ).momentum().x2());
  COMPARE(Pbox.data(box_light ).momentum().x3(), Pdef.data(def_light ).momentum().x3());
  COMPARE(Pbox.data(box_slow  ).momentum().x0(), Pdef.data(def_slow  ).momentum().x0());
  COMPARE(Pbox.data(box_slow  ).momentum().x1(), Pdef.data(def_slow  ).momentum().x1());
  COMPARE(Pbox.data(box_slow  ).momentum().x2(), Pdef.data(def_slow  ).momentum().x2());
  COMPARE(Pbox.data(box_slow  ).momentum().x3(), Pdef.data(def_slow  ).momentum().x3());
  COMPARE(Pbox.data(box_static).position().x1(), Pdef.data(def_static).position().x1());
  COMPARE(Pbox.data(box_static).position().x2(), Pdef.data(def_static).position().x2());
  COMPARE(Pbox.data(box_static).position().x3(), Pdef.data(def_static).position().x3());
  COMPARE(Pbox.data(box_light ).position().x1(), Pdef.data(def_light ).position().x1());
  COMPARE(Pbox.data(box_light ).position().x2(), Pdef.data(def_light ).position().x2());
  COMPARE(Pbox.data(box_light ).position().x3(), Pdef.data(def_light ).position().x3());
  COMPARE(Pbox.data(box_slow  ).position().x1(), Pdef.data(def_slow  ).position().x1());
  COMPARE(Pbox.data(box_slow  ).position().x2(), Pdef.data(def_slow  ).position().x2());
  COMPARE(Pbox.data(box_slow  ).position().x3(), Pdef.data(def_slow  ).position().x3());
}
