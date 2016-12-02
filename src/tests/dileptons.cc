/*
 *    Copyright (c) 2016
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "unittest.h"
#include "setup.h"

#include "../include/decayactiondilepton.h"

using namespace Smash;

TEST(init_particle_types) {
  // enable debugging output
  create_all_loggers(Configuration(""));

  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PDG\n"
      "π  0.138    7.7e-9 111 211\n"
      "e⁻ 0.000511 0      11\n"
      "γ  0        0      22\n");
}

TEST(init_decay_modes) {
  DecayModes::load_decaymodes(
      "π\n"
      "0.98823 0 γ γ\n"
      "0.01174 0  e⁻ e⁺ γ\n");
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
  const float piz_width = total_weight<DecayBranch>(type_piz.get_partial_widths(srts));
  FUZZY_COMPARE(piz_width, 7.7e-9f);
  DecayBranchPtr &mode = dil_modes[0];
  // π⁰ decay action
  const auto act = make_unique<DecayActionDilepton>(piz, 0.f,
                                                    mode->weight()/piz_width);
  act->add_decay(std::move(mode));

  // sample the final state and sum up all weights
  const int N_samples = 1E7;
  const int st = 1E6;
  float weight_sum = 0.;
  printf("sampling Dalitz ...\n");
  for (int i = 0; i < N_samples; i++) {
    if (i%st == 0) {
      std::cout << "progress (" << i/st << "/" << N_samples/st << ")" << std::endl;
    }
    act->generate_final_state();
    weight_sum += act->raw_weight_value();
  }
  std::cout << std::endl;
  std::cout << "weight_sum / N_samples = " << weight_sum / N_samples << std::endl;
  std::cout << "for # samples: " << N_samples << std::endl;
  std::cout << "should be --> 0.01174" << std::endl;
  std::cout << std::endl;

  // verify that the shining weight for the π⁰ Dalitz decay is correct
  // (to an accuracy of five percent)
  COMPARE_RELATIVE_ERROR(weight_sum / N_samples, 0.01174f, 0.05f);
}
