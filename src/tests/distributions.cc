/*
 *
 *    Copyright (c) 2014-2015,2017-2018,2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/distributions.h"

using namespace smash;

static double test_breit_wigner(const double srts, const double m,
                                const double w) {
  const double s = srts * srts;
  const double sw2 = s * w * w;
  const double smm = s - m * m;
  return 2. * s * w / (M_PI * (smm * smm + sw2));
}

TEST(breitwigner) {
  // tests the Breit-Wigner implementation for values between 0.05 and 10 GeV
  for (int s_i = 1; s_i < 200; ++s_i) {
    double s = s_i * 0.05;
    for (int m_i = 1; m_i < 200; ++m_i) {
      double m = m_i * 0.05;
      for (int w_i = 1; w_i < 200; ++w_i) {
        double w = w_i * 0.05;
        COMPARE_ABSOLUTE_ERROR(breit_wigner(s, m, w),
                               test_breit_wigner(s, m, w), 1e-6)
            << s << "/" << m << "/" << w;
      }
    }
  }
}

TEST(cauchy) {
  const double m0 = 0.770;
  const double gamma = 0.150;
  const double peak_value = 1. / (M_PI * gamma);
  // cauchy: half maximum at full width
  FUZZY_COMPARE(cauchy(m0, m0, gamma), peak_value);
  COMPARE_ABSOLUTE_ERROR(cauchy(m0 + gamma, m0, gamma), peak_value / 2., 1e-6);
  COMPARE_ABSOLUTE_ERROR(cauchy(m0 - gamma, m0, gamma), peak_value / 2., 1e-6);
  // breit_wigner_nonrel: half maximum at half width
  FUZZY_COMPARE(breit_wigner_nonrel(m0, m0, gamma), peak_value * 2.);
  COMPARE_ABSOLUTE_ERROR(breit_wigner_nonrel(m0 + gamma / 2., m0, gamma),
                         peak_value, 1e-6);
  COMPARE_ABSOLUTE_ERROR(breit_wigner_nonrel(m0 - gamma / 2., m0, gamma),
                         peak_value, 1e-6);
}
