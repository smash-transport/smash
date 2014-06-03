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
#include "../include/boxmodus.h"
#include "../include/collidermodus.h"
#include "../include/nucleusmodus.h"
#include "../include/experiment.h"
#include "../include/configuration.h"

#include <boost/filesystem.hpp>

using namespace Smash;

TEST(initialize_box) {
  Configuration conf(TEST_CONFIG_PATH);
  conf["Modi"]["Box"]["INITIAL_CONDITION"] = 1;
  conf["Modi"]["Box"]["LENGTH"] = 5.0;
  conf["Modi"]["Box"]["TEMPERATURE"] = 0.5;
  conf["Modi"]["Box"]["START_TIME"] = 0.2;
  ExperimentParameters param{{0.f, 1.f}, 1.f, 0.0, 1};
  BoxModus b(conf["Modi"], param);

  // the mass and temperature have been fixed to give a particle number
  // which is as close to an integer as possible (while still > 1).
  // The value is 131.000022, meaning this has a change of 0.0022 % of
  // giving 132 particles.
  Particles P{"smashon 0.7834 -1.0 222\n", ""};
  // should return START_TIME and set P:
  COMPARE(b.initial_conditions(&P, param), 0.2f);
  COMPARE(P.size(), 131);
  COMPARE(P.data(67).pdgcode(), 0x222);
  // we might also look at other properties of Particles, like total
  // momentum and such.
}

TEST(initialize_collider) {
  Configuration conf(TEST_CONFIG_PATH);
  conf["Modi"]["Collider"]["SQRTS"] = 1.6;
  conf["Modi"]["Collider"]["PROJECTILE"] = "222";
  conf["Modi"]["Collider"]["TARGET"] = "222";
  ExperimentParameters param{{0.f, 1.f}, 1.f, 0.0, 1};
  ColliderModus c(conf["Modi"], param);
  Particles P{"smashon 0.7834 -1.0 222\n", ""};
  COMPARE(c.initial_conditions(&P, param), -1.f);
  COMPARE(P.size(), 2);
  COMPARE(P.data(0).pdgcode(), 0x222);
  COMPARE(P.data(1).pdgcode(), 0x222);
  // here, also, we might check other properties.
}

TEST(initialize_nucleus_normal) {
  Configuration conf(TEST_CONFIG_PATH);
  conf["Modi"]["Nucleus"]["SQRTSNN"] = 1.6;
  conf.take({"Modi", "Nucleus", "Projectile"});
  conf.take({"Modi", "Nucleus", "Target"});
  conf["Modi"]["Nucleus"]["Projectile"]["PARTICLES"]["222"] = 1;
  conf["Modi"]["Nucleus"]["Target"]["PARTICLES"]["222"] = 8;
  conf["Modi"]["Nucleus"]["SQRTS_N"][0] = "222";
  conf["Modi"]["Nucleus"]["SQRTS_N"][1] = "222";
  conf["Modi"]["Nucleus"]["INITIAL_DISTANCE"] = 0;
  ExperimentParameters param{{0.f, 1.f}, 1.f, 0.0, 1};
  NucleusModus n(conf["Modi"], param);
  Particles P{"smashon 0.7834 -1.0 222", ""};
  COMPARE(n.initial_conditions(&P, param), 0.f);
  COMPARE(P.size(), 9);
  // for this, a test might include different initializations that
  // should throw an exception. This might warrant a separate file
  // nucleusmodus.cc.
}

TEST_CATCH(initialize_nucleus_low_energy, NucleusModus::InvalidEnergy) {
  Configuration conf(TEST_CONFIG_PATH);
  conf["Modi"]["Nucleus"]["SQRTSNN"] = 1.0;
  conf.take({"Modi", "Nucleus", "Projectile"});
  conf.take({"Modi", "Nucleus", "Target"});
  conf["Modi"]["Nucleus"]["Projectile"]["PARTICLES"]["222"] = 1;
  conf["Modi"]["Nucleus"]["Target"]["PARTICLES"]["222"] = 8;
  conf["Modi"]["Nucleus"]["SQRTS_N"][0] = "222";
  conf["Modi"]["Nucleus"]["SQRTS_N"][1] = "222";
  conf["Modi"]["Nucleus"]["INITIAL_DISTANCE"] = 0;
  ExperimentParameters param{{0.f, 1.f}, 1.f, 0.0, 1};
  NucleusModus n(conf["Modi"], param);
  Particles P{"smashon 0.7834 -1.0 222", ""};
  n.initial_conditions(&P, param);
}

TEST_CATCH(initialize_nucleus_empty_projectile, NucleusModus::NucleusEmpty) {
  Configuration conf(TEST_CONFIG_PATH);
  conf["Modi"]["Nucleus"]["SQRTSNN"] = 1.6;
  conf.take({"Modi", "Nucleus", "Projectile"});
  conf.take({"Modi", "Nucleus", "Target"});
  conf["Modi"]["Nucleus"]["Projectile"]["PARTICLES"]["222"] = 0;
  conf["Modi"]["Nucleus"]["Target"]["PARTICLES"]["222"] = 8;
  conf["Modi"]["Nucleus"]["SQRTS_N"][0] = 0;
  conf["Modi"]["Nucleus"]["SQRTS_N"][1] = 0;
  conf["Modi"]["Nucleus"]["INITIAL_DISTANCE"] = 0;
  ExperimentParameters param{{0.f, 1.f}, 1.f, 0.0, 1};
  NucleusModus n(conf["Modi"], param);
  Particles P{"smashon 0.7834 -1.0 222", ""};
  n.initial_conditions(&P, param);
}

TEST_CATCH(initialize_nucleus_empty_target, NucleusModus::NucleusEmpty) {
  Configuration conf(TEST_CONFIG_PATH);
  conf["Modi"]["Nucleus"]["SQRTSNN"] = 1.6;
  conf.take({"Modi", "Nucleus", "Projectile"});
  conf.take({"Modi", "Nucleus", "Target"});
  conf["Modi"]["Nucleus"]["Projectile"]["PARTICLES"]["222"] = 8;
  conf["Modi"]["Nucleus"]["Target"]["PARTICLES"]["222"] = 0;
  conf["Modi"]["Nucleus"]["SQRTS_N"][0] = 0;
  conf["Modi"]["Nucleus"]["SQRTS_N"][1] = 0;
  conf["Modi"]["Nucleus"]["INITIAL_DISTANCE"] = 0;
  ExperimentParameters param{{0.f, 1.f}, 1.f, 0.0, 1};
  NucleusModus n(conf["Modi"], param);
  Particles P{"smashon 0.7834 -1.0 222", ""};
  n.initial_conditions(&P, param);
}

// TEST(initialize_sphere) {
//   // really don't know what to do here.
// }
