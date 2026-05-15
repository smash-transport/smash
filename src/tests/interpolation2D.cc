/*
 *
 *    Copyright (c) 2020,2022,2026
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/interpolation2D.h"

#include <vector>

#include "setup.h"

using namespace smash;

static auto set_up_x_y_z_values() {
  const std::vector<double> x = {1, 2, 3, 4, 5};
  const std::vector<double> y = {1, 4, 8, 12};
  const std::vector<double> z = {1, 3, 0, 5, 0, 7, 3, 8, 9, 1,
                                 2, 5, 4, 5, 6, 1, 4, 7, 9, 2};
  return std::make_tuple(x, y, z);
}

TEST_CATCH(fail_N_points, std::invalid_argument) {
  const std::vector<double> x = {1, 2, 3, 4, 5};
  const std::vector<double> y = {0, 1, 2};
  const std::vector<double> z = {1, 2, 0, 0, 0, 0, 0, 8, 9, 1, 2, 3, 4, 5, 2};

  /* Try creating 2D interpolation with too few points, which is expected
   * to raise an exception */
  InterpolateData2DSpline(x, y, z);
}

TEST_CATCH(fail_dimensions, std::invalid_argument) {
  const std::vector<double> x = {1, 2, 3, 4, 5};
  const std::vector<double> y = {0, 1, 2, 3};
  const std::vector<double> z = {1, 2, 0, 0, 0, 8, 9, 1, 2, 3, 4, 5, 0};

  /* Try creating 2D interpolation with not-fitting dimensions, which is
   * expected to raise an exception */
  InterpolateData2DSpline(x, y, z);
}

TEST_CATCH(x_vector_not_sorted, std::invalid_argument) {
  auto [x, y, z] = set_up_x_y_z_values();
  std::swap(x[0], x[1]);

  InterpolateData2DSpline(x, y, z);
}

TEST_CATCH(y_vector_not_sorted, std::invalid_argument) {
  auto [x, y, z] = set_up_x_y_z_values();
  std::swap(y[0], y[1]);

  InterpolateData2DSpline(x, y, z);
}

TEST(interpolate_bicubic) {
  const double accuracy = 1e-6;
  const auto [x, y, z] = set_up_x_y_z_values();
  const InterpolateData2DSpline f = InterpolateData2DSpline(x, y, z);

  // check exact values at the nodes
  FUZZY_COMPARE(f(2, 4), 3.0);
  FUZZY_COMPARE(f(4, 1), 5.0);
  FUZZY_COMPARE(f(5, 8), 6.0);

  // check interpolation in x direction
  COMPARE_RELATIVE_ERROR(f(1.3, 1), 2.170375, accuracy);
  COMPARE_RELATIVE_ERROR(f(2.5, 4), 4.977678, accuracy);
  COMPARE_RELATIVE_ERROR(f(4.8, 12), 3.849142, accuracy);

  // check interpolation in y direction
  COMPARE_RELATIVE_ERROR(f(2, 6), 4.043269, accuracy);
  COMPARE_RELATIVE_ERROR(f(3, 8.1), 3.929960, accuracy);
  COMPARE_RELATIVE_ERROR(f(5, 9.5), 5.530273, accuracy);

  // check interpolation in both directions
  COMPARE_RELATIVE_ERROR(f(1.5, 5.3), 4.246599, accuracy);
  COMPARE_RELATIVE_ERROR(f(2.95, 9.1), 3.825197, accuracy);
  COMPARE_RELATIVE_ERROR(f(4.33, 11.4), 7.024691, accuracy);
}

TEST_CATCH(x_value_out_of_lower_bound, std::out_of_range) {
  const auto [x, y, z] = set_up_x_y_z_values();
  const InterpolateData2DSpline f =
      InterpolateData2DSpline(x, y, z, ExtrapolationType::None);
  f(0.5, 4);
}

TEST_CATCH(x_value_out_of_upper_bound, std::out_of_range) {
  const auto [x, y, z] = set_up_x_y_z_values();
  const InterpolateData2DSpline f =
      InterpolateData2DSpline(x, y, z, ExtrapolationType::None);
  f(5.5, 4);
}

TEST_CATCH(y_value_out_of_lower_bound, std::out_of_range) {
  const auto [x, y, z] = set_up_x_y_z_values();
  const InterpolateData2DSpline f =
      InterpolateData2DSpline(x, y, z, ExtrapolationType::None);
  f(2, 0.8);
}

TEST_CATCH(y_value_out_of_upper_bound, std::out_of_range) {
  const auto [x, y, z] = set_up_x_y_z_values();
  const InterpolateData2DSpline f =
      InterpolateData2DSpline(x, y, z, ExtrapolationType::None);
  f(2, 13);
}

TEST(extrapolate_zero) {
  // check constant extrapolation if x or y values are out of bounds
  const auto [x, y, z] = set_up_x_y_z_values();
  const InterpolateData2DSpline f =
      InterpolateData2DSpline(x, y, z, ExtrapolationType::Zero);

  // x out of bounds
  FUZZY_COMPARE(f(0.5, 4), 0.);
  FUZZY_COMPARE(f(8, 4), 0.);

  // y out of bounds
  FUZZY_COMPARE(f(2, 0.8), 0.);
  FUZZY_COMPARE(f(5, 16), 0.);

  // x and y out of bounds
  FUZZY_COMPARE(f(0.5, 0.8), 0.);
  FUZZY_COMPARE(f(7, 16), 0.);
}

TEST(extrapolate_constant) {
  // check constant extrapolation if x or y values are out of bounds
  const auto [x, y, z] = set_up_x_y_z_values();
  const InterpolateData2DSpline f =
      InterpolateData2DSpline(x, y, z, ExtrapolationType::Constant);

  // x out of bounds
  FUZZY_COMPARE(f(0.5, 4), f(1, 4));
  FUZZY_COMPARE(f(8, 4), f(5, 4));

  // y out of bounds
  FUZZY_COMPARE(f(2, 0.8), f(2, 1));
  FUZZY_COMPARE(f(5, 16), f(5, 12));

  // x and y out of bounds
  FUZZY_COMPARE(f(0.5, 0.8), f(1, 1));
  FUZZY_COMPARE(f(7, 16), f(5, 12));
}

TEST_CATCH(invalid_extrapolation_type, std::invalid_argument) {
  const auto [x, y, z] = set_up_x_y_z_values();
  InterpolateData2DSpline(x, y, z, ExtrapolationType::Linear);
}
