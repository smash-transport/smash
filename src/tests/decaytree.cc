/*
 *
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "setup.h"

#include <typeinfo>

#include "../scatteractionsfinder.cc"

using namespace smash;

TEST(init_particle_types) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PARITY PDG\n"
      "π  0.138  7.7e-9  -  111    211\n"
      "K  0.494  0       -  311    321\n"
      "ρ  0.776  0.149   -  113    213\n"
      "f₂ 1.275  0.185   +  225");
}

TEST(init_decay_channels) {
  const std::string decays_input(
      "ρ\n"
      "1.\t1\tπ π\n \n"
      "f₂\n"
      "0.50\t0\tπ π\n"
      "0.25\t0\tρ ρ\n"
      "0.20\t0\tρ ρ ρ\n"
      "0.05\t0\tK K̅\n");
  DecayModes::load_decaymodes(decays_input);
  ParticleType::check_consistency();
}

TEST(decaytree_correctness) {
  const ParticleTypePtr a = &ParticleType::find(0x211);
  const ParticleTypePtr b = &ParticleType::find(-0x211);
  const ParticleTypePtr pi0 = &ParticleType::find(0x111);
  const ParticleTypePtr rho0 = &ParticleType::find(0x113);
  const ParticleTypePtr f2 = &ParticleType::find(0x225);

  const double total_cross_section = 100.0;  // = 30 + 45 + 25
  decaytree::Node tree(a->name() + b->name(), total_cross_section, {a, b},
                       {a, b}, {a, b}, {});
  ParticleTypePtrList initial_particles1 = {a, b}, initial_particles2 = {a, b},
                      initial_particles3 = {a, b}, final_particles1 = {rho0},
                      final_particles2 = {f2}, final_particles3 = {a, b};
  auto& process_node1 =
      tree.add_action(rho0->name(), 30.0, std::move(initial_particles1),
                      std::move(final_particles1));
  decaytree::add_decays(process_node1);

  auto& process_node2 =
      tree.add_action(f2->name(), 45.0, std::move(initial_particles2),
                      std::move(final_particles2));
  decaytree::add_decays(process_node2);

  auto& process_node3 = tree.add_action(a->name() + b->name(), 25.0,
                                        std::move(initial_particles3),
                                        std::move(final_particles3));
  decaytree::add_decays(process_node3);

  tree.print();
  auto final_state_xs = tree.final_state_cross_sections();
  deduplicate(final_state_xs);

  UnitTest::setFuzzyness<double>(4);
  double xs_partial_sum = 0.0;
  for (const auto& p : final_state_xs) {
    if (p.name_ == a->name() + b->name()) {
      FUZZY_COMPARE(p.cross_section_, 70.0);
    } else if (p.name_ == pi0->name() + pi0->name()) {
      FUZZY_COMPARE(p.cross_section_, 7.5);
    } else if (p.name_ == "K̅⁻K⁺" || p.name_ == "K̅⁰K⁰") {
      FUZZY_COMPARE(p.cross_section_, 45.0 * 0.025);
    } else {
      std::cout << p.name_ << " " << p.cross_section_ << std::endl;
    }
    xs_partial_sum += p.cross_section_;
  }
  FUZZY_COMPARE(xs_partial_sum, total_cross_section);
}
