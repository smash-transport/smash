/*
 *    Copyright (c) 2016-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "unittest.h"  // This include has to be first

#include "setup.h"

#include "../include/smash/decayactiondilepton.h"

using namespace smash;

TEST(init_particle_types) {
  // enable debugging output
  create_all_loggers(Configuration(""));

  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PARITY PDG\n"
      "π  0.138  7.7e-9 - 111 211\n"
      "η  0.548 1.31e-6 - 221\n"
      "e⁻ 0.000511 0    +  11\n"
      "γ  0        0    +  22\n");
}

TEST(init_decay_modes) {
  DecayModes::load_decaymodes(
      "π\n"
      "0.98823 0 γ γ\n"
      "0.01174 0  e⁻ e⁺ γ\n"
      "\n"
      "η\n"
      "0.393   0  γ γ\n"
      "0.326   0  π⁰ π⁰ π⁰\n"
      "0.227   0  π⁺ π⁻ π⁰\n"
      "0.046   0  π⁺ π⁻ γ\n"
      "6.9e-3  0  e⁻ e⁺ γ\n");
}

TEST(pion_decay) {
  // set up a π⁰ at rest
  const ParticleType &type_piz = ParticleType::find(0x111);
  ParticleData piz{type_piz};
  piz.set_4momentum(type_piz.mass(),           // pole mass
                    ThreeVector(0., 0., 0.));  // at rest
  const auto srts = piz.effective_mass();

  // Dalitz decay π⁰ -> e⁺ e⁻ γ
  DecayBranchList dil_modes = type_piz.get_partial_widths_dilepton(srts);
  COMPARE(dil_modes.size(), 1u);
  const double piz_width =
      total_weight<DecayBranch>(type_piz.get_partial_widths(srts));
  FUZZY_COMPARE(piz_width, 7.7e-9);
  DecayBranchPtr &mode = dil_modes[0];
  // π⁰ decay action
  const auto act =
      make_unique<DecayActionDilepton>(piz, 0., mode->weight() / piz_width);
  act->add_decay(std::move(mode));

  // sample the final state and sum up all weights
  const int N_samples = 1E5;
  double weight_sum = 0.;
  printf("sampling pion Dalitz ...\n");
  for (int i = 0; i < N_samples; i++) {
    act->generate_final_state();
    weight_sum += act->get_total_weight();
  }
  // verify that the shining weight for the π⁰ Dalitz decay is correct
  // (to an accuracy of five percent)
  // the result for the π⁰ Dalitz will never match the BR exactly (even with
  // more event), because the analytic result of the formula that we use for the
  // differntial width (s. decaytype.cc) results in 4,7% overshoot of the BR
  COMPARE_RELATIVE_ERROR(weight_sum / N_samples, 0.01174, 0.05);
}

TEST(eta_decay) {
  // set up a η at rest
  const ParticleType &type_etaz = ParticleType::find(0x221);
  ParticleData etaz{type_etaz};
  etaz.set_4momentum(type_etaz.mass(),          // pole mass
                     ThreeVector(0., 0., 0.));  // at rest
  const auto srts = etaz.effective_mass();

  // Dalitz decay η -> e⁺ e⁻ γ
  DecayBranchList dil_modes = type_etaz.get_partial_widths_dilepton(srts);
  COMPARE(dil_modes.size(), 1u);
  const double etaz_width =
      total_weight<DecayBranch>(type_etaz.get_partial_widths(srts));
  FUZZY_COMPARE(etaz_width, 1.31e-6);
  DecayBranchPtr &mode = dil_modes[0];
  // π⁰ decay action
  const auto act =
      make_unique<DecayActionDilepton>(etaz, 0., mode->weight() / etaz_width);
  act->add_decay(std::move(mode));

  // sample the final state and sum up all weights
  const int N_samples = 1E5;
  double weight_sum = 0.;
  printf("sampling eta Dalitz ...\n");
  for (int i = 0; i < N_samples; i++) {
    act->generate_final_state();
    weight_sum += act->get_total_weight();
  }
  // verify that the shining weight for the η Dalitz decay is correct
  // (to an accuracy of five percent)
  COMPARE_RELATIVE_ERROR(weight_sum / N_samples, 0.0069, 0.05);
}
