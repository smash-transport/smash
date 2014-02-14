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

TEST(test_take_multiple) {
  Configuration conf(TEST_CONFIG_PATH);
  double d = conf.take({"Box", "LENGTH"});
  COMPARE(d, 10.);
  d = conf.take({"Box", "TEMPERATURE"});
  COMPARE(d, 0.2);
  int i = conf.take({"Box", "INITIAL_CONDITION"});
  COMPARE(i, 1);
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

TEST(has_value) {
  Configuration conf(TEST_CONFIG_PATH);
  VERIFY(conf.has_value({"Sphere", "RADIUS"}));
  VERIFY(conf.has_value({"Sphere", "RADIUS"}));
}

TEST(take_removes_entry) {
  Configuration conf(TEST_CONFIG_PATH);
  VERIFY(conf.has_value({"Sphere", "RADIUS"}));
  conf.take({"Sphere", "RADIUS"});
  VERIFY(!conf.has_value({"Sphere", "RADIUS"}));
}

TEST(check_unused_report) {
  std::string reference;
  Configuration conf(TEST_CONFIG_PATH);
  conf.take({"General", "MODUS"});
  conf.take({"General", "EPS"});
  conf.take({"General", "STEPS"});
  conf.take({"General", "UPDATE"});
  conf.take({"General", "RANDOMSEED"});
  conf.take({"General", "SIGMA"});
  conf.take({"General", "TESTPARTICLES"});
  conf.take({"General", "NEVENTS"});
  conf.take({"Box", "LENGTH"});
  conf.take({"Box", "TEMPERATURE"});
  conf.take({"Box", "INITIAL_CONDITION"});
  reference =
      "Collider:\n  SQRTS: 1.0\n  TARGET: -211\n  PROJECTILE: 211\nSphere:\n  "
      "RADIUS: 5.0";
  COMPARE(conf.unused_values_report(), reference);

  conf.take({"Sphere", "RADIUS"});
  reference = "Collider:\n  SQRTS: 1.0\n  TARGET: -211\n  PROJECTILE: 211";
  COMPARE(conf.unused_values_report(), reference);

  conf.take({"Collider", "PROJECTILE"});
  reference = "Collider:\n  SQRTS: 1.0\n  TARGET: -211";
  COMPARE(conf.unused_values_report(), reference);

  conf.take({"Collider", "SQRTS"});
  reference = "Collider:\n  TARGET: -211";
  COMPARE(conf.unused_values_report(), reference);

  conf.take({"Collider", "TARGET"});
  reference = "{}";
  COMPARE(conf.unused_values_report(), reference);
}
