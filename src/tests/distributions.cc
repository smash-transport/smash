/*
 *
 *    Copyright (c) 2014-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "../include/smash/distributions.h"

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

TEST(maxwell) {
  // tests the Maxwell-Boltzmann implementation for energies between 0
  // and 10 GeV and Temperatures between 0.01 and 1 GeV.
  UnitTest::setFuzzyness<double>(2);
  for (int e_i = 0; e_i < 200; ++e_i) {
    const double energy = e_i * 0.05;
    const double fourpie2 = (4.0 * M_PI) * (energy * energy);
    for (int t_i = 10; t_i < 1000; ++t_i) {
      const double temperature = t_i * 0.001;
      const double ratio = energy / temperature;
      // Avoid underflows in the exponential
      // (the estimate of maxratio does not work without the -1).
      const double maxratio = std::log(std::numeric_limits<double>::max()) - 1;
      if (ratio > maxratio) {
        continue;
      };
      COMPARE_ABSOLUTE_ERROR(
          density_integrand(energy, energy * energy, temperature),
          fourpie2 * exp(-ratio), 1e-16)
          << "E = " << energy << ", T = " << temperature;
    }
  }
}
