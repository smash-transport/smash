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
#include "../include/decayaction.h"
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
  // Create particle H, A1 and make sure their properties are as intended
  ParticleData H{ParticleType::find(0x50661)};
  ParticleData A1{ParticleType::find(0x10661)};
  const float m0_H = H.type().mass();
  const float G0_H = H.type().width_at_pole();
  const float m0_A1 = A1.type().mass();
  H.set_4momentum(H.type().mass() + 3.0f, ThreeVector(1.0, 0.0, 0.0));
  const float m_H = H.effective_mass();
  COMPARE(m0_H, 8.0f);
  COMPARE(G0_H, 1.0f);
  COMPARE(m_H, 11.0f);
  COMPARE(m0_A1, 1.0f);

  // Initialize decays of H and check their properties
  DecayBranchList H_decays = H.type().get_partial_widths(m_H);
  VERIFY(H_decays.size() == 3);
  float tmp1, tmp2, width_expected;

  int decaymodes_counter = 0;
  for (const auto &mode : H_decays) {
    float ang_mom = mode->type().angular_momentum();
    float width = mode->type().width(m0_H, G0_H, m_H);
    std::cout << "Decaymode " << decaymodes_counter << ": " <<
                 typeid(mode->type()).name() << ", " <<
                 "angular momentum: " << ang_mom <<
                 ", width: " << width << std::endl;
    VERIFY(ang_mom == 0);
    // Check if mass-dependent width behaves as expected
    switch (decaymodes_counter) {
      // Semistable two-body decay H -> A3 + A2
      case 0:
        // Here the following mathematica code was used for check:
        /* mA2 = 1.5; m0A3 = 3.0; G0A3 = 0.2; m0H = 8.0; mH = 11.0;
           mA3min = 2.5; G0H = 1.0; L = 1.6;
           PostCutoff[m_] := (L^4 + ((mA3min + mA2)^2 - m0H^2)^2/4)/
                      (L^4 + (m^2 - ((mA3min + mA2)^2 + m0H^2)/2)^2)
          rho[m_] :=  NIntegrate[
             2 mA3*( mA3*G0A3/Pi)/((mA3^2 - m0A3^2)^2 + mA3^2 G0A3^2)*
             Sqrt[(m^2 + mA2^2 - mA3^2)^2/m^2/4 - mA2^2], {mA3, mA3min,
             m - mA2} ]
          rho[mH]/rho[m0H]*PostCutoff[mH]^2*G0H  // results in 0.0122167
        */
        /* It might seem weird that resulting width is so unphysically
           small. This is because of the Post form-factor, which is intended
           for resonanse masses smaller than 2 GeV. For higher masses it
           does not give physically reasonable results. But this is only
           code test, so we can live with it.
        */
        COMPARE_RELATIVE_ERROR(width, 0.0122167f, 1.e-7);
        break;
      // Stable two-body decay H -> A1 + A1
      case 1:
        /* In the case of equal product masses expressions simplify
           and one gets Gamma(m) / Gamma(m0) =
           = \sqrt((1 - 4 m_{A1}^2 / m_H^2) / (1 - 4 m_{A1}^2 / m0_H^2))
        */
        tmp1 = 2 * m0_A1 / m_H;
        tmp2 = 2 * m0_A1 / m0_H;
        width_expected = G0_H * std::sqrt((1.f - tmp1*tmp1)/(1.f - tmp2*tmp2));
        COMPARE_RELATIVE_ERROR(width, width_expected, 1.e-7);
        break;
      // three-body decay H -> A2 + A2 + A1
      case 2:
        COMPARE_RELATIVE_ERROR(width, G0_H, 1.e-7);
        break;
      // Should never get here
      default: VERIFY(0 == 1);
    }
    decaymodes_counter++;
  }
  VERIFY(decaymodes_counter == 3);

  const float time_of_execution = 4.5f;
  const auto act = make_unique<DecayAction>(H, time_of_execution);
  std::cout << *act << std::endl;
}
