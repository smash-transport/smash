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
#include "../include/spheremodus.h"

#include <boost/filesystem.hpp>

using namespace Smash;

TEST(init_particle_types) {
  ParticleType::create_type_list(
      "smashon 0.4 0.0 661\n");
}

TEST(initialize_box) {
  Configuration conf(TEST_CONFIG_PATH);
  conf["Modi"]["Box"]["INITIAL_CONDITION"] = 1;
  conf["Modi"]["Box"]["LENGTH"] = 7.9615;
  conf["Modi"]["Box"]["TEMPERATURE"] = 0.5;
  conf["Modi"]["Box"]["START_TIME"] = 0.2;
  conf.take({"Modi", "Box", "INIT_MULTIPLICITIES"});
  conf["Modi"]["Box"]["INIT_MULTIPLICITIES"]["661"] = 724;
  ExperimentParameters param{{0.f, 1.f}, 1.f, 0.0, 1};
  BoxModus b(conf["Modi"], param);

  Particles P;
  // should return START_TIME and set P:
  COMPARE(b.initial_conditions(&P, param), 0.2f);
  COMPARE(P.size(), 724);
  COMPARE(P.data(67).pdgcode(), 0x661);
  // we might also look at other properties of Particles, like total
  // momentum and such.
  FourVector momentum(0.0, 0.0, 0.0, 0.0);
  for (auto p : P.data()) {
    momentum += p.momentum();
    VERIFY(p.position().x1() <  7.9615);
    VERIFY(p.position().x1() >= 0.0);
    VERIFY(p.position().x2() <  7.9615);
    VERIFY(p.position().x2() >= 0.0);
    VERIFY(p.position().x3() <  7.9615);
    VERIFY(p.position().x3() >= 0.0);
  }
  COMPARE_ABSOLUTE_ERROR(momentum.x1(), 0.0, 1e-12);
  COMPARE_ABSOLUTE_ERROR(momentum.x2(), 0.0, 1e-12);
  COMPARE_ABSOLUTE_ERROR(momentum.x3(), 0.0, 1e-12);
}

TEST(initialize_collider) {
  Configuration conf(TEST_CONFIG_PATH);
  conf["Modi"]["Collider"]["SQRTS"] = 1.6;
  conf["Modi"]["Collider"]["PROJECTILE"] = "661";
  conf["Modi"]["Collider"]["TARGET"] = "661";
  ExperimentParameters param{{0.f, 1.f}, 1.f, 0.0, 1};
  ColliderModus c(conf["Modi"], param);
  Particles P;
  COMPARE(c.initial_conditions(&P, param), -1.f);
  COMPARE(P.size(), 2);
  COMPARE(P.data(0).pdgcode(), 0x661);
  COMPARE(P.data(1).pdgcode(), 0x661);
  VERIFY (P.data(0).position().x1() > 0.0);
  COMPARE(P.data(1).position().x1(), 0.0);
  COMPARE(P.data(0).position().x2(), 0.0);
  COMPARE(P.data(1).position().x2(), 0.0);
  COMPARE(P.data(0).position().x3(), -P.data(1).position().x3());
  VERIFY (P.data(0).position().x3() < 0.0);
  // total momentum is checked in two steps:
  FourVector momentum_total = P.data(0).momentum() + P.data(1).momentum();
  // energy component with a possible deviation,
  COMPARE_RELATIVE_ERROR(momentum_total.x0(), 1.6, 1e-7);
  // and the spatial components exactly.
  COMPARE(momentum_total, FourVector(momentum_total.x0(), 0.0, 0.0, 0.0));
}

TEST_CATCH(initialize_collider_low_energy, ModusDefault::InvalidEnergy) {
  Configuration conf(TEST_CONFIG_PATH);
  conf["Modi"]["Collider"]["SQRTS"] = 0.5;
  conf["Modi"]["Collider"]["PROJECTILE"] = "661";
  conf["Modi"]["Collider"]["TARGET"] = "661";
  ExperimentParameters param{{0.f, 1.f}, 1.f, 0.0, 1};
  ColliderModus c(conf["Modi"], param);
}

TEST(initialize_nucleus_normal) {
  Configuration conf(TEST_CONFIG_PATH);
  conf["Modi"]["Nucleus"]["SQRTSNN"] = 1.6;
  conf.take({"Modi", "Nucleus", "Projectile"});
  conf.take({"Modi", "Nucleus", "Target"});
  conf["Modi"]["Nucleus"]["Projectile"]["PARTICLES"]["661"] = 1;
  conf["Modi"]["Nucleus"]["Target"]["PARTICLES"]["661"] = 8;
  conf["Modi"]["Nucleus"]["SQRTS_REPS"][0] = "661";
  conf["Modi"]["Nucleus"]["SQRTS_REPS"][1] = "661";
  conf["Modi"]["Nucleus"]["INITIAL_DISTANCE"] = 0;
  conf["Modi"]["Nucleus"]["Impact"]["VALUE"] = 0;
  ExperimentParameters param{{0.f, 1.f}, 1.f, 0.0, 1};
  NucleusModus n(conf["Modi"], param);
  Particles P;
  COMPARE(n.initial_conditions(&P, param), 0.f);
  COMPARE(P.size(), 9);
  for (auto p : P.data()) {
    // velocity should be +- sqrt(3/4)
    COMPARE_RELATIVE_ERROR(p.velocity().sqr(), 0.75, 1e-6);
    // this is the mass squared
    COMPARE_RELATIVE_ERROR(p.momentum().sqr(), 0.16, 1e-6);
    COMPARE(p.position().x0(), 0.0);
    COMPARE(p.pdgcode(), PdgCode(0x661));
    COMPARE_RELATIVE_ERROR(p.momentum().x0(), 0.8, 1e-6);
    COMPARE_ABSOLUTE_ERROR(p.momentum().x1(), 0.0, 1e-6);
    COMPARE_ABSOLUTE_ERROR(p.momentum().x2(), 0.0, 1e-6);
    COMPARE_RELATIVE_ERROR(std::abs(p.momentum().x3()), std::sqrt(0.48), 1e-6);
  }
  // all other things can only be tested with statistics.
}

TEST_CATCH(initialize_nucleus_low_energy, ModusDefault::InvalidEnergy) {
  Configuration conf(TEST_CONFIG_PATH);
  conf["Modi"]["Nucleus"]["SQRTSNN"] = 0.5;
  conf.take({"Modi", "Nucleus", "Projectile"});
  conf.take({"Modi", "Nucleus", "Target"});
  conf["Modi"]["Nucleus"]["Projectile"]["PARTICLES"]["661"] = 1;
  conf["Modi"]["Nucleus"]["Target"]["PARTICLES"]["661"] = 8;
  conf["Modi"]["Nucleus"]["SQRTS_REPS"][0] = "661";
  conf["Modi"]["Nucleus"]["SQRTS_REPS"][1] = "661";
  conf["Modi"]["Nucleus"]["INITIAL_DISTANCE"] = 0;
  ExperimentParameters param{{0.f, 1.f}, 1.f, 0.0, 1};
  NucleusModus n(conf["Modi"], param);
  Particles P;
  n.initial_conditions(&P, param);
}

TEST_CATCH(initialize_nucleus_empty_projectile, NucleusModus::NucleusEmpty) {
  Configuration conf(TEST_CONFIG_PATH);
  conf["Modi"]["Nucleus"]["SQRTSNN"] = 1.6;
  conf.take({"Modi", "Nucleus", "Projectile"});
  conf.take({"Modi", "Nucleus", "Target"});
  conf["Modi"]["Nucleus"]["Projectile"]["PARTICLES"]["661"] = 0;
  conf["Modi"]["Nucleus"]["Target"]["PARTICLES"]["661"] = 8;
  conf["Modi"]["Nucleus"]["SQRTS_REPS"][0] = 0;
  conf["Modi"]["Nucleus"]["SQRTS_REPS"][1] = 0;
  conf["Modi"]["Nucleus"]["INITIAL_DISTANCE"] = 0;
  ExperimentParameters param{{0.f, 1.f}, 1.f, 0.0, 1};
  NucleusModus n(conf["Modi"], param);
  Particles P;
  n.initial_conditions(&P, param);
}

TEST_CATCH(initialize_nucleus_empty_target, NucleusModus::NucleusEmpty) {
  Configuration conf(TEST_CONFIG_PATH);
  conf["Modi"]["Nucleus"]["SQRTSNN"] = 1.6;
  conf.take({"Modi", "Nucleus", "Projectile"});
  conf.take({"Modi", "Nucleus", "Target"});
  conf["Modi"]["Nucleus"]["Projectile"]["PARTICLES"]["661"] = 8;
  conf["Modi"]["Nucleus"]["Target"]["PARTICLES"]["661"] = 0;
  conf["Modi"]["Nucleus"]["SQRTS_REPS"][0] = 0;
  conf["Modi"]["Nucleus"]["SQRTS_REPS"][1] = 0;
  conf["Modi"]["Nucleus"]["INITIAL_DISTANCE"] = 0;
  ExperimentParameters param{{0.f, 1.f}, 1.f, 0.0, 1};
  NucleusModus n(conf["Modi"], param);
  Particles P;
  n.initial_conditions(&P, param);
}

TEST(initialize_sphere) {
 Configuration conf(TEST_CONFIG_PATH);
  conf["Modi"]["Sphere"]["RADIUS"] = 10;
  conf["Modi"]["Sphere"]["NUMBEROFPARTICLES"] = 300;
  conf["Modi"]["Sphere"]["SPHERETEMPERATURE"] = 0.5;
  conf["Modi"]["Sphere"]["START_TIME"] = 0.0;
  ExperimentParameters param{{0.f, 1.f}, 1.f, 0.0, 1};
  SphereModus s(conf["Modi"], param);
  Particles P;
//Is the correct number of particles in the map?
  COMPARE(s.initial_conditions(&P, param), 0.0f);
  COMPARE(P.size(), 300);
  COMPARE(P.data(67).pdgcode(), 0x661);
// total momentum check
  FourVector momentum(0.0, 0.0, 0.0, 0.0);
// position less than radius?
  float radius = 0.0;
  for (auto p : P.data()) {
    momentum += p.momentum();
    radius = sqrt(p.position().x1()*p.position().x1()+
    p.position().x2()*p.position().x2()+p.position().x3()*p.position().x3()); 
    VERIFY(radius <  10.0);
  }
  COMPARE_ABSOLUTE_ERROR(momentum.x1(), 0.0, 1e-12);
  COMPARE_ABSOLUTE_ERROR(momentum.x2(), 0.0, 1e-12);
  COMPARE_ABSOLUTE_ERROR(momentum.x3(), 0.0, 1e-12);
}
