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

TEST(init_particle_types) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PDG\n"
      "pi0 0.1350 -1.0 111\n"
      "pi+ 0.1396 -1.0 211\n"
      "rho0 0.7755 0.149 113\n"
      "rho+ 0.7755 0.149 213\n"
      "eta 0.5479 1.0e-6 221\n"
      "omega 0.7827 0.0085 223\n"
      "p 0.9383 -1.0 2212\n"
      "n 0.9396 -1.0 2112\n"
      "Delta++ 1.232 0.117 2224\n"
      "Delta+ 1.232 0.117 2214\n"
      "Delta0 1.232 0.117 2114\n"
      "Delta- 1.232 0.117 1114\n");
}

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
