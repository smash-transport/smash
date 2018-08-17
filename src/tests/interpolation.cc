/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "../include/smash/interpolation.h"

#include <vector>

using namespace smash;

TEST(interpolate_linear) {
  const auto f = InterpolateLinear<double>(0, 0, 1, 2);
  COMPARE(f(0), 0);
  COMPARE(f(1), 2);
  COMPARE(f(0.5), 1);
}

TEST(permutation) {
  const std::vector<double> x = {0, 7, 5, 4, 8, 6, 2, 1, 3, 9};
  const std::vector<double> y = {9, 2, 4, 5, 1, 3, 7, 8, 6, 0};
  const auto p =
      generate_sort_permutation(x, [&](double a, double b) { return a < b; });
  const std::vector<double> sorted_x = apply_permutation(x, p);
  const std::vector<double> correctly_sorted_x = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  COMPARE(sorted_x, correctly_sorted_x);
  const std::vector<double> permuted_y = apply_permutation(y, p);
  const std::vector<double> correctly_permuted_y = {9, 8, 7, 6, 5,
                                                    4, 3, 2, 1, 0};
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
  const std::vector<double> data = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
  COMPARE(find_index(data, -1.0), 0ul);
  COMPARE(find_index(data, 0.2), 0ul);
  COMPARE(find_index(data, 0.3), 1ul);
  COMPARE(find_index(data, 0.4), 1ul);
  COMPARE(find_index(data, 0.5), 2ul);
  COMPARE(find_index(data, 10.0), 5ul);
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
