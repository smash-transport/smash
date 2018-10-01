/*
 *
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "../include/smash/constants.h"
#include "../include/smash/kinematics.h"

using namespace smash;

TEST(center_of_velocity_v) {
  const double s = 2.9 * 2.9;
  const double ma = 0.6;
  const double mb = 1.2;
  const double v = center_of_velocity_v(s, ma, mb);
  const double gamma = 1.0 / std::sqrt(1.0 - v * v);
  const double E = gamma * (ma + mb);
  const double p = gamma * v * (ma - mb);
  FUZZY_COMPARE(s, E * E - p * p);
}

TEST(fixed_target_projectile_v) {
  const double s = 2.9 * 2.9;
  const double ma = 0.6;
  const double mb = 1.2;
  const double v = fixed_target_projectile_v(s, ma, mb);
  const double gamma = 1.0 / std::sqrt(1.0 - v * v);
  const double E = gamma * ma + mb;
  const double p = gamma * v * ma;
  UnitTest::setFuzzyness<double>(2);
  FUZZY_COMPARE(s, E * E - p * p);
}

TEST(pCM) {
  const double srts = 2.9;
  const double ma = 0.6;
  const double mb = 1.2;
  const double pcm_sqr = pCM_sqr(srts, ma, mb);
  const double pcm = pCM(srts, ma, mb);
  FUZZY_COMPARE(std::sqrt(pcm_sqr + ma * ma) + std::sqrt(pcm_sqr + mb * mb),
                srts);
  COMPARE(pcm * pcm, pcm_sqr);
  COMPARE(pCM_sqr_from_s(srts * srts, ma, mb), pcm_sqr);
  COMPARE(pCM_from_s(srts * srts, ma, mb), pcm);
}

TEST(plab_from_s_NN) {
  const double s = 2.9 * 2.9;
  const double mN = 0.938;
  const double plab = plab_from_s(s);
  const double E = std::sqrt(plab * plab + mN * mN) + mN;
  // Requiring more precision makes test fail
  COMPARE_RELATIVE_ERROR(s, E * E - plab * plab, 1.e-7);
  COMPARE_RELATIVE_ERROR(s, s_from_plab(plab, mN, mN), 1.e-7);
}

TEST(plab_from_s_KN) {
  const double s = 2.9 * 2.9;
  const double plab = plab_from_s(s, kaon_mass, nucleon_mass);
  const double E =
      std::sqrt(plab * plab + kaon_mass * kaon_mass) + nucleon_mass;
  COMPARE_RELATIVE_ERROR(s, E * E - plab * plab, 1.e-7);
  COMPARE_RELATIVE_ERROR(s, s_from_plab(plab, kaon_mass, nucleon_mass), 1.e-7);
}

TEST(plab_from_s_KN_small) {
  // At this value plab should vanish, but the function is very steep there.
  const double s = (kaon_mass + nucleon_mass) * (kaon_mass + nucleon_mass);
  // We add a small constant to avoid numerical issue with the assert.
  COMPARE_ABSOLUTE_ERROR(plab_from_s(s + 1e-9, kaon_mass, nucleon_mass), 0.0,
                         3e-4);
  // std::cout << plab_from_s(1.7, kaon_mass, nucleon_mass) << std::endl;
}

TEST(s_from_Ekin) {
  const double Ekin = 1.1;
  const double ma = 0.6;
  const double mb = 1.2;
  const double E = Ekin + ma + mb;
  const double p_sqr = (Ekin + ma) * (Ekin + ma) - ma * ma;
  const double s = s_from_Ekin(Ekin, ma, mb);
  FUZZY_COMPARE(s, E * E - p_sqr);
}

TEST(s_from_plab) {
  const double plab = 1.1;
  const double ma = 0.6;
  const double mb = 1.2;
  const double E = std::sqrt(plab * plab + ma * ma) + mb;
  UnitTest::setFuzzyness<double>(2);
  FUZZY_COMPARE(s_from_plab(plab, ma, mb), E * E - plab * plab);
}
