/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include "../include/integrate.h"

TEST(no_arguments) {
  Smash::Integrator integrate;
  for (int i = 0; i < 10; ++i) {
    const auto result = integrate(0, i, [](double) { return 1.; });
    COMPARE_ABSOLUTE_ERROR(result.value(), double(i), result.error());
  }
  for (int i = 0; i < 10; ++i) {
    const auto result = integrate(0, i, [](double x) { return x; });
    COMPARE_ABSOLUTE_ERROR(result.value(), i * i * 0.5, result.error());
  }
}

TEST(with_lambda_captures) {
  Smash::Integrator integrate;
  for (int i = 0; i < 10; ++i) {
    const auto result = integrate(0, i, [i](double x) { return x + i; });
    COMPARE_ABSOLUTE_ERROR(result.value(), i * i * 1.5, result.error()) << "i = " << i;
  }
  for (int i = 0; i < 10; ++i) {
    double y = i * 2.;
    const auto result = integrate(0, i, [&](double x) { return x * y + i; });
    COMPARE_ABSOLUTE_ERROR(result.value(), i * i + i * i * y * 0.5, result.error()) << "\ni = " << i;
  }
}

TEST(result_conversion_to_value) {
  Smash::Integrator integrate;
  for (int i = 0; i < 10; ++i) {
    const auto result1 = integrate(0, i, [](double x) { return x; });
    const double result2 = integrate(0, i, [](double x) { return x; });
    COMPARE(result1.value(), result2);
  }
}
