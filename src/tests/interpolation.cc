/*
 *
 *    Copyright (c) 2015-2018,2020,2022,2026
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/interpolation.h"

#include <vector>

using namespace smash;

static InterpolateDataLinear<double> set_up_interpolate_data_linear_sorted(
    ExtrapolationType extrapolation_type) {
  const std::vector<double> x = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  const std::vector<double> y = {1, 2, 0, 0, 0, 0, 0, 8, 9};
  return InterpolateDataLinear<double>(x, y, extrapolation_type);
}

static InterpolateDataSpline set_up_interpolate_data_spline_sorted(
    ExtrapolationType extrapolation_type) {
  const std::vector<double> x = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  const std::vector<double> y = x;
  return InterpolateDataSpline(x, y, extrapolation_type);
}

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

TEST(find_index) {
  const std::vector<double> data = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
  COMPARE(find_index(data, -1.0), 0ul);
  COMPARE(find_index(data, 0.2), 0ul);
  COMPARE(find_index(data, 0.3), 1ul);
  COMPARE(find_index(data, 0.4), 1ul);
  COMPARE(find_index(data, 0.5), 2ul);
  COMPARE(find_index(data, 10.0), 5ul);
}

TEST(interpolate_data_linear) {
  const InterpolateDataLinear<double> f =
      set_up_interpolate_data_linear_sorted(ExtrapolationType::None);
  COMPARE(f(1.5), 1.5);
  COMPARE(f(5), 0.0);
}

TEST_CATCH(interpolate_data_linear_value_out_of_lower_bound,
           std::out_of_range) {
  const InterpolateDataLinear<double> f =
      set_up_interpolate_data_linear_sorted(ExtrapolationType::None);
  f(0);
}

TEST_CATCH(interpolate_data_linear_value_out_of_upper_bound,
           std::out_of_range) {
  const InterpolateDataLinear<double> f =
      set_up_interpolate_data_linear_sorted(ExtrapolationType::None);
  f(10);
}

TEST(interpolate_data_linear_unsorted) {
  std::vector<double> x = {5, 6, 7, 9, 3, 4, 1, 2, 8};
  std::vector<double> y = {0, 0, 0, 9, 0, 0, 1, 2, 8};
  const InterpolateDataLinear<double> f(x, y);
  x.clear();
  y.clear();
  COMPARE(f(1.5), 1.5);
  COMPARE(f(5), 0.0);
}

TEST(interpolate_data_linear_with_zero_extrapolation) {
  const InterpolateDataLinear<double> f =
      set_up_interpolate_data_linear_sorted(ExtrapolationType::Zero);
  COMPARE(f(1.5), 1.5);
  COMPARE(f(0), 0.);
  COMPARE(f(10), 0.);
}

TEST(interpolate_data_linear_with_constant_extrapolation) {
  const InterpolateDataLinear<double> f =
      set_up_interpolate_data_linear_sorted(ExtrapolationType::Constant);
  COMPARE(f(1.5), 1.5);
  COMPARE(f(0), 1.);
  COMPARE(f(10), 9.);
}

TEST(interpolate_data_linear_with_linear_extrapolation) {
  const InterpolateDataLinear<double> f =
      set_up_interpolate_data_linear_sorted(ExtrapolationType::Linear);
  COMPARE(f(1.5), 1.5);
  COMPARE(f(0), 0.);
  COMPARE(f(10), 10.);
}

TEST(interpolate_data_spline) {
  const InterpolateDataSpline f =
      set_up_interpolate_data_spline_sorted(ExtrapolationType::None);
  COMPARE(f(1.5), 1.5);
}

TEST_CATCH(interpolate_data_spline_value_out_of_lower_bound,
           std::out_of_range) {
  const InterpolateDataSpline f =
      set_up_interpolate_data_spline_sorted(ExtrapolationType::None);
  f(0);
}

TEST_CATCH(interpolate_data_spline_value_out_of_upper_bound,
           std::out_of_range) {
  const InterpolateDataSpline f =
      set_up_interpolate_data_spline_sorted(ExtrapolationType::None);
  f(10);
}

TEST(interpolate_data_spline_unsorted) {
  std::vector<double> x = {7, 5, 4, 8, 6, 2, 1, 3, 9};
  std::vector<double> y = x;
  const InterpolateDataSpline f(x, y);
  x.clear();
  y.clear();
  COMPARE(f(1.5), 1.5);
}

TEST(interpolate_data_spline_with_zero_extrapolation) {
  const InterpolateDataSpline f =
      set_up_interpolate_data_spline_sorted(ExtrapolationType::Zero);
  COMPARE(f(1.5), 1.5);
  COMPARE(f(0), 0.);
  COMPARE(f(10), 0.);
}

TEST(interpolate_data_spline_with_constant_extrapolation) {
  const InterpolateDataSpline f =
      set_up_interpolate_data_spline_sorted(ExtrapolationType::Constant);
  COMPARE(f(1.5), 1.5);
  COMPARE(f(0), 1.);
  COMPARE(f(10), 9.);
}

TEST_CATCH(interpolate_data_spline_invalid_extrapolation_type,
           std::invalid_argument) {
  set_up_interpolate_data_spline_sorted(ExtrapolationType::Linear);
}
