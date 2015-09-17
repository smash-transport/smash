/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include "unittest.h"
#include "../include/interpolation.h"

#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

TEST(interpolate_linear) {
  const auto f = InterpolateLinear<double>(0, 0, 1, 2);
  COMPARE(f(0), 0);
  COMPARE(f(1), 2);
  COMPARE(f(0.5), 1);
}

TEST(permutation) {
  const std::vector<double> x = {0, 7, 5, 4, 8, 6, 2, 1, 3, 9};
  const std::vector<double> y = {9, 2, 4, 5, 1, 3, 7, 8, 6, 0};
  const auto p = generate_sort_permutation(
      x, [&](double a, double b) {
        return a < b;
      });
  const std::vector<double> sorted_x = std::move(apply_permutation(x, p));
  const std::vector<double> correctly_sorted_x = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  COMPARE(sorted_x, correctly_sorted_x);
  const std::vector<double> permuted_y = std::move(apply_permutation(y, p));
  const std::vector<double> correctly_permuted_y = {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
  COMPARE(permuted_y, correctly_permuted_y);
}

TEST(interpolate_data_linear) {
  std::vector<double> x = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  std::vector<double> y = {1, 2, 0, 0, 0, 0, 0, 8, 9};
  InterpolateDataLinear<double> f(x, y);
  x.resize(0);
  y.resize(0);
  COMPARE(f(1.5), 1.5);
  COMPARE(f(5), 0.0);
  COMPARE(f(0), 0.0);
  COMPARE(f(10), 10.0);
}

TEST(interpolate_data_linear_unsorted) {
  std::vector<double> x = {5, 6, 7, 9, 3, 4, 1, 2, 8};
  std::vector<double> y = {0, 0, 0, 9, 0, 0, 1, 2, 8};
  InterpolateDataLinear<double> f(x, y);
  x.resize(0);
  y.resize(0);
  COMPARE(f(1.5), 1.5);
  COMPARE(f(5), 0.0);
  COMPARE(f(0), 0.0);
  COMPARE(f(10), 10.0);
}

TEST(find_index) {
  const std::vector<float> data = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 };
  COMPARE(find_index(data, -1.0f), 0ul);
  COMPARE(find_index(data, 0.2f), 0ul);
  COMPARE(find_index(data, 0.3f), 1ul);
  COMPARE(find_index(data, 0.4f), 1ul);
  COMPARE(find_index(data, 0.5f), 2ul);
  COMPARE(find_index(data, 10.0f), 5ul);
}

TEST(interpolate_data_spline) {
  std::vector<double> x = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  std::vector<double> y = x;
  const InterpolateDataSpline f(x, y);
  x.resize(0);
  y.resize(0);
  COMPARE(f(1.5), 1.5);
  COMPARE(f(0), 1.0);
  COMPARE(f(10), 9.0);
}

TEST(interpolate_data_spline_unsorted) {
  std::vector<double> x = {7, 5, 4, 8, 6, 2, 1, 3, 9};
  std::vector<double> y = x;
  const InterpolateDataSpline f(x, y);
  x.resize(0);
  y.resize(0);
  COMPARE(f(1.5), 1.5);
  COMPARE(f(0), 1.0);
  COMPARE(f(10), 9.0);
}

TEST(interpolate_data_bspline) {
  std::vector<double> x = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  std::vector<double> y = x;
  const InterpolateDataBSpline f(x, y, 5);
  x.resize(0);
  y.resize(0);
  UnitTest::setFuzzyness<double>(8);
  FUZZY_COMPARE(f(1.5), 1.5);
  //FUZZY_COMPARE(f(0), 1.0);
  //FUZZY_COMPARE(f(10), 9.0);
}

TEST(interpolate_data_bspline_example) {
  constexpr size_t n = 200;
  constexpr size_t ncoeffs = 12;
  constexpr size_t k = 4;
  constexpr size_t nbreak = ncoeffs + 2 - k;
  constexpr double x_min = 0;
  constexpr double x_max = 15;

  // Generate data with noise.
  std::vector<double> x;
  x.reserve(n);
  std::vector<double> y;
  y.reserve(n);
  std::vector<double> w;
  w.reserve(n);
  gsl_rng_env_setup();
  auto r = gsl_rng_alloc(gsl_rng_default);
  for (size_t i = 0; i < n; i++) {
    const double xi = x_min + (x_max / (n - 1)) * i;
    double yi = cos(xi) * exp(-0.1 * xi);
    const double sigma = 0.1 * yi;
    const double dy = gsl_ran_gaussian(r, sigma);
    yi += dy;
    x.push_back(xi);
    y.push_back(yi);
    w.push_back(1.0 / (sigma * sigma));
  }

  // Fit.
  const auto f = InterpolateDataBSpline(x, y, nbreak);

  // Output smoothed curve.
  for (double xi = x_min; xi < x_max; xi += 0.1) {
    printf("%f %f\n", xi, f(xi));
  }

  gsl_rng_free(r);
}
