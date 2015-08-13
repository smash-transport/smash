/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include "unittest.h"

#include "../include/kinematics.h"

using namespace Smash;

TEST(center_of_velocity_v) {
  const float s = 2.9f*2.9f;
  const float ma = 0.6f;
  const float mb = 1.2f;
  const double v = center_of_velocity_v(s, ma, mb);
  const double gamma = 1.0 / std::sqrt(1.0 - v*v);
  const float E = gamma * (ma + mb);
  const float p = gamma * v * (ma - mb);
  FUZZY_COMPARE(s, E*E -p*p);
}

TEST(fixed_target_projectile_v) {
  const float s = 2.9f*2.9f;
  const float ma = 0.6f;
  const float mb = 1.2f;
  const double v = fixed_target_projectile_v(s, ma, mb);
  const double gamma = 1.0 / std::sqrt(1.0 - v*v);
  const float E = gamma*ma + mb;
  const float p = gamma*v*ma;
  FUZZY_COMPARE(s, E*E - p*p);
}

TEST(pCM) {
  const float srts = 2.9f;
  const float ma = 0.6f;
  const float mb = 1.2f;
  const float pcm_sqr = pCM_sqr(srts, ma, mb);
  const float pcm = pCM(srts, ma, mb);
  FUZZY_COMPARE(std::sqrt(pcm_sqr + ma*ma) +
                std::sqrt(pcm_sqr + mb*mb), srts);
  COMPARE(pcm*pcm, pcm_sqr);
  COMPARE(pCM_sqr_from_s(srts*srts, ma, mb), pcm_sqr);
  COMPARE(pCM_from_s(srts*srts, ma, mb), pcm);
}

TEST(plab_from_sNN) {
  const double s = 2.9*2.9;
  const double mN = 0.938;
  const double plab = plab_from_s_NN(s);
  const double E = std::sqrt(plab*plab + mN*mN) + mN;
  // Requiring more precision makes test fail
  COMPARE_RELATIVE_ERROR(s, E*E - plab*plab, 1.e-7);
  COMPARE_RELATIVE_ERROR(s, s_from_plab(plab, mN, mN), 1.e-7);
}

TEST(s_from_Ekin) {
  const double Ekin = 1.1;
  const double ma = 0.6f;
  const double mb = 1.2f;
  const double E = Ekin + ma + mb;
  const double p_sqr = (Ekin + ma) * (Ekin + ma) - ma*ma;
  const double s = s_from_Ekin(Ekin, ma, mb);
  FUZZY_COMPARE(s, E*E - p_sqr);
}

TEST(s_from_plab) {
  const double plab = 1.1;
  const double ma = 0.6f;
  const double mb = 1.2f;
  const double E = std::sqrt(plab*plab + ma*ma) + mb;
  FUZZY_COMPARE(s_from_plab(plab, ma, mb), E*E - plab*plab);
}
