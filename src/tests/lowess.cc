/*
 *
 *    Copyright (c) 2015-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "../include/smash/lowess.h"

#include <numeric>

TEST(smooth_trivial) {
  std::vector<double> x(20);
  std::iota(x.begin(), x.end(), 0);
  std::vector<double> y(20);
  for (size_t i = 0; i < 10; i++) {
    y[i] = 0;
  }
  for (size_t i = 10; i < 20; i++) {
    y[i] = 1;
  }
  const auto result = smash::smooth(x, y, 0.4, 3, 0.);
  for (size_t i = 0; i < 20; i++) {
    COMPARE_ABSOLUTE_ERROR(result[i], y[i], 1e-7);
  }
}

TEST(smooth_flat) {
  std::vector<double> x(20);
  std::iota(x.begin(), x.end(), 0);
  std::vector<double> y(20);
  for (size_t i = 0; i < 20; i++) {
    y[i] = 0;
  }
  const auto result = smash::smooth(x, y, 0.4, 3, 0.);
  for (size_t i = 0; i < 20; i++) {
    COMPARE_ABSOLUTE_ERROR(result[i], 0., 1e-7);
  }
}

TEST(smooth_range) {
  std::vector<double> x(20);
  std::iota(x.begin(), x.end(), 0);
  std::vector<double> y = x;
  // This test fails for 3 robustness iterations!
  constexpr size_t robustness_iterations = 0;
  const auto result = smash::smooth(x, y, 0.4, robustness_iterations, 0.);
  for (size_t i = 0; i < 20; i++) {
    COMPARE_ABSOLUTE_ERROR(result[i], x[i], 1e-7);
  }
}

TEST(smooth_simple) {
  std::vector<double> x = {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                           10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
  std::vector<double> y = {-0.76741118, 0.69245631,  2.39950921,  2.53647578,
                           2.32918222,  5.6595567,   6.66367639,  4.95611415,
                           8.8123281,   10.45977518, 11.21428038, 12.29296866,
                           12.78028477, 12.7597147,  13.78278698, 15.24549405,
                           16.25987014, 16.09290966, 16.54311784, 18.68219495};
  std::vector<double> smoothed_y = {
      -0.626034455349546, 0.56507171201094, 1.75962718897954, 2.95796332584499,
      4.15606361537761,   5.34733969366442, 6.52229821799894, 7.70815938803622,
      8.87590555190343,   9.940975860297,   10.8981138457948, 11.7851424727769,
      12.6188717296918,   13.409849737403,  14.1516996584552, 14.9180658146586,
      15.695660019874,    16.4783034134255, 17.2617441530539, 18.0459201716397};
  const auto result = smash::smooth(x, y, 2. / 3, 3, 0.);
  for (size_t i = 0; i < 20; i++) {
    COMPARE_ABSOLUTE_ERROR(result[i], smoothed_y[i], 1e-7);
  }
}
