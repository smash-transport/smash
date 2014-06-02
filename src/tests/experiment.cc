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
