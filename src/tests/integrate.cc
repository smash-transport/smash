/*
 *
 *    Copyright (c) 2015,2017-2018,2020,2022,2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/integrate.h"

/* Test one-dimensional integration: using simple integration function, the
 * result is expected to be few ULP away from the analytical result. We tolerate
 * 3 at most here. The comparison with the absolute error calculated by GSL
 * cannot be done "at one sigma" and therefore we double the error in the
 * absolute comparison. Note that the Integrator class has an accuracy_absolute_
 * member which by default is way larger then the GSL returned error in these
 * simple cases. */
constexpr int fuzzyness = 3;

TEST(one_dim_no_arguments) {
  // This test performse used algorithm sometimes underestimates the true error
  // by a few bits of precision.
  vir::test::setFuzzyness<double>(fuzzyness);
  smash::Integrator integrate;
  for (int i = 0; i < 10; ++i) {
    const auto result = integrate(0, i, [](double) { return 1.; });
    FUZZY_COMPARE(result.value(), double(i));
    COMPARE_ABSOLUTE_ERROR(result.value(), double(i), 2 * result.error());
  }
  for (int i = 0; i < 10; ++i) {
    const auto result = integrate(0, i, [](double x) { return x; });
    FUZZY_COMPARE(result.value(), i * i * 0.5);
    COMPARE_ABSOLUTE_ERROR(result.value(), i * i * 0.5, 2 * result.error());
  }
}

TEST(one_dim_with_lambda_captures) {
  vir::test::setFuzzyness<double>(fuzzyness);
  smash::Integrator integrate;
  for (int i = 0; i < 10; ++i) {
    const auto result = integrate(0, i, [i](double x) { return x + i; });
    FUZZY_COMPARE(result.value(), i * i * 1.5);
    COMPARE_ABSOLUTE_ERROR(result.value(), i * i * 1.5, 2 * result.error())
        << "i = " << i;
  }
  for (int i = 0; i < 10; ++i) {
    const double y = i * 2.;
    const auto result = integrate(0, i, [&](double x) { return x * y + i; });
    const double exact_result = i * i * (i + 1);
    FUZZY_COMPARE(result.value(), exact_result);
    COMPARE_ABSOLUTE_ERROR(result.value(), exact_result, 2 * result.error())
        << "\ni = " << i;
  }
}

TEST(one_dim_result_conversion_to_value) {
  smash::Integrator integrate;
  for (int i = 0; i < 10; ++i) {
    const auto result1 = integrate(0, i, [](double x) { return x; });
    const double result2 = integrate(0, i, [](double x) { return x; });
    COMPARE(result1.value(), result2);
  }
}

// test two-dimensional Monte-Carlo integration
TEST(two_dim) {
  smash::Integrator2d integrate;
  /* Here the errors are statistical. We check that the values agree within
   * 4 sigma. It is very unlikely that they lie outside the 4 sigma band. */
  constexpr int Nsigma = 4;
  // constant integrand
  for (int i = 0; i < 10; ++i) {
    const auto result =
        integrate(0, i, 0, i, [](double, double) { return 1.; });
    COMPARE_ABSOLUTE_ERROR(result.value(), double(i * i),
                           Nsigma * result.error());
  }
  // linear only in one dim
  for (int i = 0; i < 10; ++i) {
    const auto result =
        integrate(0, i, 0, i, [](double x, double) { return x; });
    COMPARE_ABSOLUTE_ERROR(result.value(), i * i * i * 0.5,
                           Nsigma * result.error());
  }
  // linear in both dims (factorizable)
  for (int i = 0; i < 10; ++i) {
    const auto result =
        integrate(0, i, 0, i, [](double x, double y) { return x * y; });
    COMPARE_ABSOLUTE_ERROR(result.value(), i * i * i * i * 0.25,
                           Nsigma * result.error());
  }
  // non-factorizable
  for (int i = 0; i < 10; ++i) {
    const auto result = integrate(
        0, i, 0, i, [](double x, double y) { return std::sqrt(x + y); });
    COMPARE_ABSOLUTE_ERROR(
        result.value(), 8. / 15. * (2. * std::sqrt(2.) - 1.) * pow(i, 5. / 2.),
        Nsigma * result.error() * 5);
  }
}
