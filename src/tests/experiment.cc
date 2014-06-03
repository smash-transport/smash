/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"

#include "../include/experiment.h"
#include "../include/configuration.h"

#include <boost/filesystem.hpp>

using namespace Smash;

TEST(create_box) {
  Configuration conf(TEST_CONFIG_PATH);
  conf["General"]["MODUS"] = "Box";
  VERIFY(!!ExperimentBase::create(conf));
}

TEST(create_collider) {
  Configuration conf(TEST_CONFIG_PATH);
  conf["General"]["MODUS"] = "Collider";
  VERIFY(!!ExperimentBase::create(conf));
}

TEST(create_nucleus) {
  Configuration conf(TEST_CONFIG_PATH);
  conf["General"]["MODUS"] = "Nucleus";
  VERIFY(!!ExperimentBase::create(conf));
}

// TEST(create_sphere) {
//   Configuration conf(TEST_CONFIG_PATH);
//   conf["General"]["MODUS"] = "Sphere";
//   VERIFY(!!ExperimentBase::create(conf));
// }
//
TEST_CATCH(create_invalid, ExperimentBase::InvalidModusRequest) {
  Configuration conf(TEST_CONFIG_PATH);
  conf["General"]["MODUS"] = "Invalid";
  ExperimentBase::create(conf);
}

// TEST(experiment_parameters) {
//   Configuration conf(TEST_CONFIG_PATH);
//   conf["General"]["TESTPARTICLES"] = 2;
//   conf["General"]["SIGMA"] = 3;
//   conf["General"]["DELTA_TIME"] = 4;
//   conf["General"]["OUTPUT_INTERVAL"] = 5;
//   ExperimentParameters param = create_experiment_parameters(conf);
//   COMPARE(param.output_interval, 5.f);
//   COMPARE(param.cross_section, 3.0 / 2.0);
//   COMPARE(param.testparticles, 2);
// }
