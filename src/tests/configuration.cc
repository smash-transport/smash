/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "tests/unittest.h"
#include "include/configuration.h"

#include <boost/filesystem.hpp>

TEST(create_object) {
  Configuration conf(TEST_CONFIG_PATH);
}

TEST(check_config_general_contents) {
  Configuration conf(TEST_CONFIG_PATH);

  const auto general = conf["General"];
  COMPARE(general["MODUS"        ].as<std::string>(), "Collider");
  COMPARE(general["EPS"          ].as<double     >(), 0.01);
  COMPARE(general["STEPS"        ].as<int        >(), 1000);
  COMPARE(general["UPDATE"       ].as<int        >(), 10);
  COMPARE(general["RANDOMSEED"   ].as<int        >(), 1);
  COMPARE(general["SIGMA"        ].as<double     >(), 10.0);
  COMPARE(general["TESTPARTICLES"].as<int        >(), 1);
  COMPARE(general["NEVENTS"      ].as<int        >(), 1);
}

TEST(check_config_collider_contents) {
  Configuration conf(TEST_CONFIG_PATH);
  const auto collider = conf["Collider"];
  COMPARE(collider["PROJECTILE"].as<int        >(), 211);
  COMPARE(collider["TARGET"    ].as<int        >(), -211);
  COMPARE(collider["SQRTS"     ].as<double     >(), 1.0);
}

TEST(test_take) {
  Configuration conf(TEST_CONFIG_PATH);
  double d = conf.take({"Sphere", "RADIUS"});
  COMPARE(d, 5.);
}

TEST_CATCH(take_incorrect_type, Configuration::IncorrectTypeInAssignment) {
  Configuration conf(TEST_CONFIG_PATH);
  int i = conf.take({"Sphere", "RADIUS"});
  COMPARE(i, 5);
}

TEST(take_always_converts_to_string) {
  Configuration conf(TEST_CONFIG_PATH);
  std::string s = conf.take({"Sphere", "RADIUS"});
  COMPARE(s, "5.0");
}
