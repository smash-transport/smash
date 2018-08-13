/*
 *
 *    Copyright (c) 2014-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "../include/smash/particledata.h"
#include "../include/smash/processbranch.h"

using namespace smash;

TEST(init_particle_types) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PDG\n"
      "σ    1.1 1.1 9876542\n"
      "σino 1.1 1.1 1234568\n");
}

TEST(assign_default) {
  CollisionBranch branch(0., ProcessType::StringSoftSingleDiffractiveAX);
  FUZZY_COMPARE(branch.weight(), 0.);
  COMPARE(branch.get_type(), ProcessType::StringSoftSingleDiffractiveAX);
}
TEST(assign_1_particle) {
  PdgCode smashon("9876542");
  CollisionBranch branch(ParticleType::find(smashon), 1.234,
                         ProcessType::Elastic);
  FUZZY_COMPARE(branch.weight(), 1.234);
}
TEST(assign_2_particle) {
  PdgCode smashon("9876542");
  CollisionBranch branch(ParticleType::find(smashon),
                         ParticleType::find(smashon), 2.345,
                         ProcessType::Elastic);
  FUZZY_COMPARE(branch.weight(), 2.345);
}

TEST(lists) {
  const ParticleType &smashon(ParticleType::find(PdgCode("9876542")));
  const ParticleType &smashino(ParticleType::find(PdgCode("1234568")));
  CollisionBranch branch(smashon, smashino, 2.345, ProcessType::Elastic);
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
}

TEST(add_particle) {
  std::vector<ParticleTypePtr> list = {
      &ParticleType::find(PdgCode("9876542")),
      &ParticleType::find(PdgCode("1234568")),
      &ParticleType::find(PdgCode("-1234568")),
  };
  CollisionBranch branch(list, 1.2, ProcessType::Elastic);
  COMPARE(branch.particle_types().size(), 3u);
}

TEST(weights) {
  CollisionBranch branch(0., ProcessType::Elastic);
  branch.set_weight(0.34);
  COMPARE(branch.weight(), 0.34);
  // double is intentional here.
  branch.set_weight(0.33);
  COMPARE(branch.weight(), 0.33);
}
