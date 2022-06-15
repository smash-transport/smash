/*
 *    Copyright (c) 2020-2022
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
      "Λ  1.116   0         +      3122\n"
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
  act1 = std::make_unique<ScatterActionMulti>(incoming_pis, 0.05);
  ScatterActionMultiPtr act2;
  act2 = std::make_unique<ScatterActionMulti>(incoming_2pi_eta, 0.05);

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
  act1 = std::make_unique<ScatterActionMulti>(incoming_pi, 0.05);
  ScatterActionMultiPtr act2;
  act2 = std::make_unique<ScatterActionMulti>(incoming_n, 0.05);

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

TEST(threebody_integral_I3) {
  // Make sure incoming particles got their masses set
  // calculate_I3 uses effective mass , not type mass
  ParticleData p{ParticleType::find(0x2212)};  // p
  p.set_4momentum(p.type().mass(), ThreeVector(0.1, 0.2, 0.3));
  ParticleData n{ParticleType::find(0x2112)};  // n
  n.set_4momentum(n.type().mass(), ThreeVector(0.1, 0.2, 0.3));
  ParticleData pip{ParticleType::find(0x211)};  // pi+
  pip.set_4momentum(pip.type().mass(), ThreeVector(0.1, 0.2, 0.3));

  const ParticleList& incoming_piNN{pip, p, n};
  ScatterActionMultiPtr act_piNN =
      std::make_unique<ScatterActionMulti>(incoming_piNN, 0.0);

  std::vector<double> srts{2.5, 3.0, 4.0, 5.0, 6.0};
  // Answers from Mathematica
  std::vector<double> correct_I3{1.21877, 7.21727, 49.8706, 168.18, 414.645};
  VERIFY(srts.size() == correct_I3.size());
  for (size_t i = 0; i < srts.size(); i++) {
    COMPARE_RELATIVE_ERROR(act_piNN->calculate_I3(srts[i]), correct_I3[i],
                           1e-4);
  }
}

TEST(phi4_parametrization) {
  ParticleData N{ParticleType::find(0x2212)};   // p
  ParticleData pi{ParticleType::find(0x211)};   // pi+
  ParticleData La{ParticleType::find(0x3122)};  // Lambda

  const ParticleList& incoming_piNNN = {pi, N, N, N};
  const ParticleList& incoming_NNNN = {N, N, N, N};
  const ParticleList& incoming_piNNLa = {pi, N, N, La};
  const ParticleList& incoming_NNNLa = {N, N, N, La};

  ScatterActionMultiPtr act_piNNN =
      std::make_unique<ScatterActionMulti>(incoming_piNNN, 0.05);
  ScatterActionMultiPtr act_NNNN =
      std::make_unique<ScatterActionMulti>(incoming_NNNN, 0.05);
  ScatterActionMultiPtr act_piNNLa =
      std::make_unique<ScatterActionMulti>(incoming_piNNLa, 0.05);
  ScatterActionMultiPtr act_NNNLa =
      std::make_unique<ScatterActionMulti>(incoming_NNNLa, 0.05);

  const double srts = 4.5;  // GeV
  const double s = srts * srts;
  const double piNNN = act_piNNN->parametrizaton_phi4(s);
  const double NNNN = act_NNNN->parametrizaton_phi4(s);
  const double piNNLa = act_piNNLa->parametrizaton_phi4(s);
  const double NNNLa = act_NNNLa->parametrizaton_phi4(s);
  // Expectations are computed by numerical integration in Mathematica
  // Declared relatice precision of parametrizations is 10^-3
  COMPARE_RELATIVE_ERROR(piNNN, 3.18511e-6, 1e-3) << piNNN;
  COMPARE_RELATIVE_ERROR(NNNN, 3.62849e-7, 1e-3) << NNNN;
  COMPARE_RELATIVE_ERROR(piNNLa, 2.08782e-6, 1e-3) << piNNLa;
  COMPARE_RELATIVE_ERROR(NNNLa, 1.4423e-7, 1e-3) << NNNLa;
}
