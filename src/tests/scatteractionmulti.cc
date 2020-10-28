/*
 *    Copyright (c) 2016-2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include <vir/test.h>  // This include has to be first

#include "setup.h"

#include "../include/smash/scatteractionmulti.h"

using namespace smash;
using smash::Test::Momentum;
using smash::Test::Position;

const FourVector pos_a = Position{0., -0.1, 0., 0.};
const FourVector pos_b = Position{0., 0.1, 0., 0.};
const FourVector pos_c = Position{0., 0.0, 0., 0.};

TEST(init_particle_types) {
  // enable debugging output
  create_all_loggers(Configuration(""));

  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PARITY PDG\n"
      "π  0.138   7.7e-9    -      111      211\n"
      "N  0.938   0         +      2112     2212\n"
      "η  0.548   1.31e-6   -      221\n"
      "η' 0.958   1.96e-4   -      331\n"
      "φ  1.019   4.25e-3   -      333\n"
      "ω  0.783   8.49e-3   -      223\n"
      "d  1.8756  0         +      1000010020\n"
      "d' 1.8856  0.100     -      1000010021\n");
}

TEST(init_decay_modes) {
  DecayModes::load_decaymodes(
      "ω\n"
      "0.893    1  π⁺ π⁻ π⁰\n"
      "\n"
      "φ\n"
      "0.152   1  π⁺ π⁻ π⁰\n"
      "\n"
      "η'\n"
      "0.426  0  π⁺ π⁻ η\n"
      "0.228  0  π⁰ π⁰ η\n"
      "\n"
      "d'\n"
      "1.      1   N N\n");
  ParticleType::check_consistency();
  sha256::Hash hash;
  hash.fill(0);
  IsoParticleType::tabulate_integrals(hash, "");
}

TEST(three_meson_to_one) {
  Momentum some_momentum{1.1, 1.0, 0., 0.};

  ParticleData pip{ParticleType::find(0x211)};  // pi+
  pip.set_4momentum(some_momentum);

  ParticleData pim{ParticleType::find(-0x211)};  // pi-
  pim.set_4momentum(some_momentum);

  ParticleData piz{ParticleType::find(0x111)};  // pi0
  piz.set_4momentum(some_momentum);

  ParticleData eta{ParticleType::find(0x221)};  // eta
  eta.set_4momentum(some_momentum);

  const ParticleList& incoming_pis = {pip, pim, piz};
  const ParticleList& incoming_2pi_eta = {pip, pim, eta};

  MultiParticleReactionsBitSet incl_all_multi_set =
      MultiParticleReactionsBitSet().set();

  ScatterActionMultiPtr act1;
  act1 = make_unique<ScatterActionMulti>(incoming_pis, 0.05);
  ScatterActionMultiPtr act2;
  act2 = make_unique<ScatterActionMulti>(incoming_2pi_eta, 0.05);

  act1->add_possible_reactions(0.1, 8.0, incl_all_multi_set);
  act2->add_possible_reactions(0.1, 8.0, incl_all_multi_set);

  // Phi and omega 3->1 reactions should be found for three different pions
  COMPARE(act1->reaction_channels().size(), 2u);
  // Only eta-prime 3->1 reaction should be found for 2pi and eta
  COMPARE(act2->reaction_channels().size(), 1u);

  // Verify that process type are correct
  VERIFY(act1->reaction_channels()[0]->get_type() ==
         ProcessType::MultiParticleThreeMesonsToOne);
  VERIFY(act1->reaction_channels()[1]->get_type() ==
         ProcessType::MultiParticleThreeMesonsToOne);
  VERIFY(act2->reaction_channels()[0]->get_type() ==
         ProcessType::MultiParticleThreeMesonsToOne);
}

TEST(deuteron_three_to_two) {
  Momentum some_momentum{1.1, 1.0, 0., 0.};
  Momentum some_other_momentum{1.1, -1.0, 0., 0.};

  ParticleData p{ParticleType::find(0x2212)};  // p
  p.set_4momentum(some_momentum);

  ParticleData n{ParticleType::find(0x2112)};  // n
  n.set_4momentum(some_momentum);

  ParticleData pipp{ParticleType::find(0x211)};  // pi+
  pipp.set_4momentum(some_other_momentum);

  ParticleData antip{ParticleType::find(-0x2212)};  // anti-p
  antip.set_4momentum(some_other_momentum);

  const ParticleList& incoming_pi = {pipp, p, n};
  const ParticleList& incoming_n = {antip, p, n};

  MultiParticleReactionsBitSet incl_all_multi_set =
      MultiParticleReactionsBitSet().set();

  ScatterActionMultiPtr act1;
  act1 = make_unique<ScatterActionMulti>(incoming_pi, 0.05);
  ScatterActionMultiPtr act2;
  act2 = make_unique<ScatterActionMulti>(incoming_n, 0.05);

  act1->add_possible_reactions(0.1, 8.0, incl_all_multi_set);
  act2->add_possible_reactions(0.1, 8.0, incl_all_multi_set);

  // πpn → πd should be found
  COMPARE(act1->reaction_channels().size(), 1u);
  // N̅np → N̅d should be found
  COMPARE(act2->reaction_channels().size(), 1u);

  // Verify that process type are correct
  VERIFY(act1->reaction_channels()[0]->get_type() ==
         ProcessType::MultiParticleThreeToTwo);
  VERIFY(act2->reaction_channels()[0]->get_type() ==
         ProcessType::MultiParticleThreeToTwo);
}
