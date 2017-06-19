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
      "Λ 3.000 0.3 50661\n"
      "η1⁰ 0.400 -1.0 10661\n"
      "η2⁰ 0.600 -1.0 20661\n"
      "η3⁰ 1.200 0.2 30661");
}

TEST(init_decay_channels) {
  // A3 -> A1 + A2
  // H -> A3 + A2, A1 + A1, A2 + A2 + A1
  const std::string decays_input(
      "η3 \n"
      "2.0\t0\tη2 η1\n \n"
      "Λ\n \n"
      " 0.5 \t0\tη3 η2\n \n"
      " 1.0 \t0\tη1 η1\n \n"
      " 1.5 \t0\tη2⁰ η1⁰ η1⁰\n");
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
  H.set_4momentum(H.type().mass() + 1.0f, ThreeVector(1.0, 0.0, 0.0));
  const float m_H = H.effective_mass();
  COMPARE(m0_H, 3.0f);
  COMPARE(G0_H, 0.3f);
  COMPARE(m_H,  4.0f);
  COMPARE(m0_A1, 0.4f);
  // Check consistency for width at pole
  COMPARE(H.type().total_width(m0_H), G0_H);

  // Initialize decays of H and check their properties
  DecayBranchList H_decays = H.type().get_partial_widths(m_H);
  COMPARE(H_decays.size(), 3);
  float tmp1, tmp2, width_expected;

  int decaymodes_counter = 0;
  for (const auto &mode : H_decays) {
    const float ang_mom = mode->type().angular_momentum();
    const float width = mode->weight();
    const DecayType &type = mode->type();
    std::cout << "Decaymode " << decaymodes_counter << ": " <<
                 typeid(type).name() << ", " <<
                 "angular momentum: " << ang_mom <<
                 ", width: " << width << std::endl;
    COMPARE(ang_mom, 0);
    // Check if mass-dependent width behaves as expected
    switch (decaymodes_counter) {
      // Semistable two-body decay H -> A3 + A2
      case 0:
        /* Result obtained with MATHEMATICA code:
           mA2 = 0.9; mA1 = 0.5; m0A3 = 1.0; G0A3 = 0.2; m0H = 3.0; mH = 4.0;
           mA3min = 1.4; G0H = 0.3; L = 1.6; s0 = mA3min + mA2;
           PostCutoff[m_] := (L^4 + (s0^2 - m0H^2)^2/4)/(L^4 + (m^2 - (s0^2 + m0H^2)/2)^2)
           PfOverM[m_, m1_, m2_] := Sqrt[(m^2 + m1^2 - m2^2)^2/(4 m^2) - m1^2]/m
           GA3[mA3_] := G0A3*PfOverM[mA3, mA1, mA2]/PfOverM[m0A3, mA1, mA2]
           A[m_] := 2/Pi m^2 GA3[m]/((m^2 - m0A3^2)^2 + (m GA3[m])^2)
           rho[m_] :=  NIntegrate[A[mA3]*PfOverM[mA3, mA1, mA2], {mA3, mA3min, m - mA2}]
           rho[mH]/rho[m0H]*PostCutoff[mH]^2*G0H // returns 0.00824107*/
        /* It might seem weird that resulting width is so unphysically
           small. This is because of the Post form-factor, which is intended
           for resonanse masses smaller than 2 GeV. For higher masses it
           does not give physically reasonable results. But this is only
           code test, so we can live with it.
        */
        COMPARE_RELATIVE_ERROR(width, 0.00824107f / 6.f, 5.e-2);
        break;
      // Stable two-body decay H -> A1 + A1
      case 1:
        /* In the case of equal product masses expressions simplify
           and one gets Gamma(m) / Gamma(m0) =
           = \sqrt((1 - 4 m_{A1}^2 / m_H^2) / (1 - 4 m_{A1}^2 / m0_H^2))
        */
        tmp1 = 2 * m0_A1 / m_H;
        tmp2 = 2 * m0_A1 / m0_H;
        width_expected = G0_H / 3. * std::sqrt((1.f - tmp1*tmp1)/(1.f - tmp2*tmp2));
        COMPARE_RELATIVE_ERROR(width, width_expected, 1.e-6);
        break;
      // three-body decay H -> A2 + A2 + A1
      case 2:
        COMPARE_RELATIVE_ERROR(width, G0_H / 2.f, 1.e-7);
        break;
      // Should never get here
      default: VERIFY(false);
    }
    decaymodes_counter++;
  }
  COMPARE(decaymodes_counter, 3);

  const float time_of_execution = 4.5f;
  const auto act = make_unique<DecayAction>(H, time_of_execution);
  std::cout << *act << std::endl;
}
