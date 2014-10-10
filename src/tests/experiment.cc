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

namespace particles_txt {
#include <particles.txt.h>
}

using namespace Smash;

TEST(init_particle_types) {
  ParticleType::create_type_list(particles_txt::data);
}

TEST(create_box) {
  Configuration conf(TEST_CONFIG_PATH);
  conf["General"]["Modus"] = "Box";
  VERIFY(!!ExperimentBase::create(conf));
}

TEST(create_collider) {
  Configuration conf(TEST_CONFIG_PATH);
  conf["General"]["Modus"] = "Collider";
  VERIFY(!!ExperimentBase::create(conf));
}

TEST(create_nucleus) {
  Configuration conf(TEST_CONFIG_PATH);
  conf["General"]["Modus"] = "Nucleus";
  VERIFY(!!ExperimentBase::create(conf));
}

TEST(create_sphere) {
  Configuration conf(TEST_CONFIG_PATH);
  conf["General"]["Modus"] = "Sphere";
  VERIFY(!!ExperimentBase::create(conf));
}

TEST_CATCH(create_invalid, ExperimentBase::InvalidModusRequest) {
  Configuration conf(TEST_CONFIG_PATH);
  conf["General"]["Modus"] = "Invalid";
  ExperimentBase::create(conf);
}
