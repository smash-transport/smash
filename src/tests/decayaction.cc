/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <typeinfo>

#include "unittest.h"
#include "setup.h"

#include "../include/decaymodes.h"
#include "../include/action.h"
#include "../include/cxx14compat.h"

using namespace Smash;

TEST(init_particle_types) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PDG\n"
      "H 8.000 1.0 50661\n"
      "A1 1.000 -1.0 10661\n"
      "A2 1.500 -1.0 20661\n"
      "A3 3.000 0.2 30661");
}

TEST(init_decay_channels) {
  // A3 -> A1 + A2
  // H -> A3 + A2, A1 + A1, A2 + A2 + A1
  const std::string decays_input(
      "30661 \t# A3\n"
      "2.0\t0\t20661 10661\t# A2 A1 \n \n"
      " 50661\t# H\n \n"
      " 0.5 \t0\t30661 20661\t# A3 A2 \n \n"
      " 1.0 \t0\t10661 10661\t# A1 A1 \n \n"
      " 1.5 \t0\t20661 20661 10661\t# A2 A2 A1\n");
  DecayModes::load_decaymodes(decays_input);
  ParticleType::check_consistency();
}

TEST(create_decayaction) {
  ParticleData H{ParticleType::find(0x50661)};
  H.set_4momentum(H.type().mass() + 0.1, ThreeVector(1.0, 0.0, 0.0));
  const float time_of_execution = 4.5f;
  // H.type() - ParticleType
  // H.type().decay_modes() - DecayModes
  // H.type().decay_modes().decay_mode_list() - DecayBranchList
  int decaymodes_counter = 0;
  for (const auto &mode : H.type().decay_modes().decay_mode_list()) {
    // mode - unique pointer to DecayBranch
    // mode->type() - DecayType
    std::cout << "Decaymode " << decaymodes_counter << ": " <<
                 typeid(mode->type()).name() << ", " <<
                 "angular momentum: " << mode->type().angular_momentum() <<
                 ", width: " << mode->type().width(H.type().mass(),
                                                   H.type().width_at_pole(),
                                                    H.effective_mass())
                 << std::endl;
    decaymodes_counter++;
  }
  DecayBranchList H_decays =
          H.type().get_partial_widths(H.effective_mass());

  const auto act = make_unique<DecayAction>(H, time_of_execution);
  std::cout << *act << std::endl;
}
