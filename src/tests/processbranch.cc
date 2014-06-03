/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include "../include/processbranch.h"
#include "../include/particledata.h"

using namespace Smash;

TEST(assign_default) {
  ProcessBranch branch;
  FUZZY_COMPARE(branch.weight(), -1.f);
  FUZZY_COMPARE(branch.type(), 0);
}
TEST(assign_1_particle) {
  PdgCode smashon("9876543");
  ProcessBranch branch(smashon, 1.234f, 8);
  FUZZY_COMPARE(branch.weight(), 1.234f);
  FUZZY_COMPARE(branch.type(), 8);
}
TEST(assign_2_particle) {
  PdgCode smashon("9876543");
  ProcessBranch branch(smashon, smashon, 2.345, 9);
  FUZZY_COMPARE(branch.weight(), 2.345f);
  FUZZY_COMPARE(branch.type(), 9);
}

TEST(lists) {
  PdgCode smashon("9876543");
  PdgCode smashino("1234567");
  ProcessBranch branch(smashon, smashino, 2.345, 9);
  std::vector<PdgCode> list = branch.pdg_list();
  COMPARE(list.size(), 2);
  COMPARE(list.at(0), smashon);
  COMPARE(list.at(1), smashino);

  ParticleList particles = branch.particle_list();
  COMPARE(particles.size(), 2);
  COMPARE(particles.at(0).pdgcode(), smashon);
  COMPARE(particles.at(1).pdgcode(), smashino);
  COMPARE(particles.at(0).momentum().x0(), 0.0);
  COMPARE(particles.at(0).momentum().x1(), 0.0);
  COMPARE(particles.at(0).momentum().x2(), 0.0);
  COMPARE(particles.at(0).momentum().x3(), 0.0);
  COMPARE(particles.at(1).momentum().x0(), 0.0);
  COMPARE(particles.at(1).momentum().x1(), 0.0);
  COMPARE(particles.at(1).momentum().x2(), 0.0);
  COMPARE(particles.at(1).momentum().x3(), 0.0);

  branch.clear();

  std::vector<PdgCode> new_list = branch.pdg_list();
  COMPARE(new_list.size(), 0);
}

TEST(add_particle) {
  PdgCode smashon("9876543");
  PdgCode smashino("1234567");
  PdgCode anti_smashino("-1234567");
  ProcessBranch branch(smashon, smashino, 1.2, 1);
  branch.add_particle(anti_smashino);
  COMPARE(branch.pdg_list().size(), 3);
}
TEST(set_particles) {
  PdgCode smashon("9876543");
  PdgCode smashino("1234567");
  PdgCode anti_smashino("-1234567");
  std::vector<PdgCode> list = {smashon, smashino, anti_smashino};
  ProcessBranch branch;
  branch.set_particles(list);
  COMPARE(branch.pdg_list().size(), 3);
}

TEST(weights) {
  ProcessBranch branch;
  branch.set_weight(0.34f);
  COMPARE(branch.weight(), 0.34f);
  // double is intentional here.
  branch.set_weight(0.33);
  COMPARE(branch.weight(), 0.33f);
}

TEST(types) {
  ProcessBranch branch;
  branch.set_type(18);
  COMPARE(branch.type(), 18);
}
