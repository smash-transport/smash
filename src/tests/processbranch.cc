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

TEST(init_particle_types) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PDG\n"
      "smashon 1.1 1.1 9876542\n"
      "smashino 1.1 1.1 1234568\n"
      "anto_smashino 1.1 1.1 -1234568\n");
}

TEST(assign_default) {
  ProcessBranch branch;
  FUZZY_COMPARE(branch.weight(), -1.f);
}
TEST(assign_1_particle) {
  PdgCode smashon("9876542");
  ProcessBranch branch(ParticleType::find(smashon), 1.234f);
  FUZZY_COMPARE(branch.weight(), 1.234f);
}
TEST(assign_2_particle) {
  PdgCode smashon("9876542");
  ProcessBranch branch(ParticleType::find(smashon), ParticleType::find(smashon), 2.345);
  FUZZY_COMPARE(branch.weight(), 2.345f);
}

TEST(lists) {
  const ParticleType &smashon(ParticleType::find({"9876542"}));
  const ParticleType &smashino(ParticleType::find({"1234568"}));
  ProcessBranch branch(smashon, smashino, 2.345);
  const auto &list = branch.particle_types();
  COMPARE(list.size(), 2u);
  COMPARE(list.at(0), &smashon);
  COMPARE(list.at(1), &smashino);

  ParticleList particles = branch.particle_list();
  COMPARE(particles.size(), 2u);
  COMPARE(particles.at(0).pdgcode(), smashon.pdgcode());
  COMPARE(particles.at(1).pdgcode(), smashino.pdgcode());
  COMPARE(particles.at(0).momentum().x0(), 0.0);
  COMPARE(particles.at(0).momentum().x1(), 0.0);
  COMPARE(particles.at(0).momentum().x2(), 0.0);
  COMPARE(particles.at(0).momentum().x3(), 0.0);
  COMPARE(particles.at(1).momentum().x0(), 0.0);
  COMPARE(particles.at(1).momentum().x1(), 0.0);
  COMPARE(particles.at(1).momentum().x2(), 0.0);
  COMPARE(particles.at(1).momentum().x3(), 0.0);

  branch.clear();

  const auto &new_list = branch.particle_types();
  COMPARE(new_list.size(), 0u);
}

TEST(add_particle) {
  std::vector<ParticleTypePtr> list = {
      &ParticleType::find({"9876542"}), &ParticleType::find({"1234568"}),
      &ParticleType::find({"-1234568"}),
  };
  ProcessBranch branch(list, 1.2);
  COMPARE(branch.particle_types().size(), 3u);
}

TEST(weights) {
  ProcessBranch branch;
  branch.set_weight(0.34f);
  COMPARE(branch.weight(), 0.34f);
  // double is intentional here.
  branch.set_weight(0.33);
  COMPARE(branch.weight(), 0.33f);
}
