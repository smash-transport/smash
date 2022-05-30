/*
 *
 *    Copyright (c) 2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <vir/test.h>  // This include has to be first

#include <vector>
#include "setup.h"

#include "../include/smash/interpolation2D.h"

using namespace smash;

std::unique_ptr<InterpolateData2DSpline> interp = nullptr;
const double accuracy = 1e-6;

TEST(fail_N_points) {
  std::vector<double> x = {1, 2, 3, 4, 5};
  std::vector<double> y = {1, 2, 0};
  std::vector<double> z = {1, 2, 0, 0, 0, 0, 0, 8, 9, 1, 2, 3, 4, 5, 2};

  // Try creating 2D interpolation with too few points, which is expected
  // to raise an exception
  vir::test::expect_failure();
  interp = std::make_unique<InterpolateData2DSpline>(x, y, z);
}

TEST(fail_dimensions) {
  std::vector<double> x = {1, 2, 3, 4, 5};
  std::vector<double> y = {1, 2, 0};
  std::vector<double> z = {1, 2, 0, 0, 0, 8, 9, 1, 2, 3, 4, 5, 0};

  // Try creating 2D interpolation with not-fitting dimensions, which is
  // expected to raise an exception
  vir::test::expect_failure();
  interp = std::make_unique<InterpolateData2DSpline>(x, y, z);
}

TEST(interpolate_bicubic) {
  std::vector<double> x = {1, 2, 3, 4, 5};
  std::vector<double> y = {1, 4, 8, 12};
  std::vector<double> z = {1, 3, 0, 5, 0, 7, 3, 8, 9, 1,
                           2, 5, 4, 5, 6, 1, 4, 7, 9, 2};
  interp = std::make_unique<InterpolateData2DSpline>(x, y, z);

  // check exact values at the nodes
  FUZZY_COMPARE((*interp)(2, 4), 3.0);
  FUZZY_COMPARE((*interp)(4, 1), 5.0);
  FUZZY_COMPARE((*interp)(5, 8), 6.0);

  // check interpolation in x direction
  COMPARE_RELATIVE_ERROR((*interp)(1.3, 1), 2.170375, accuracy);
  COMPARE_RELATIVE_ERROR((*interp)(2.5, 4), 4.977678, accuracy);
  COMPARE_RELATIVE_ERROR((*interp)(4.8, 12), 3.849142, accuracy);

  // check interpolation in y direction
  COMPARE_RELATIVE_ERROR((*interp)(2, 6), 4.043269, accuracy);
  COMPARE_RELATIVE_ERROR((*interp)(3, 8.1), 3.929960, accuracy);
  COMPARE_RELATIVE_ERROR((*interp)(5, 9.5), 5.530273, accuracy);

  // check interpolation in both directions
  COMPARE_RELATIVE_ERROR((*interp)(1.5, 5.3), 4.246599, accuracy);
  COMPARE_RELATIVE_ERROR((*interp)(2.95, 9.1), 3.825197, accuracy);
  COMPARE_RELATIVE_ERROR((*interp)(4.33, 92.04), 7.553750, accuracy);
}

TEST(extrapolate_constant) {
  // check constant extrapolation if x or y values are out of bounds
  std::vector<double> x = {1, 2, 3, 4, 5};
  std::vector<double> y = {1, 4, 8, 12};
  std::vector<double> z = {1, 3, 0, 5, 0, 7, 3, 8, 9, 1,
                           2, 5, 4, 5, 6, 1, 4, 7, 9, 2};
  interp = std::make_unique<InterpolateData2DSpline>(x, y, z);

  // x out of bounds
  FUZZY_COMPARE((*interp)(0.5, 4), (*interp)(1, 4));
  FUZZY_COMPARE((*interp)(8, 4), (*interp)(5, 4));

  // y out of bounds
  FUZZY_COMPARE((*interp)(2, 0.8), (*interp)(2, 1));
  FUZZY_COMPARE((*interp)(5, 16), (*interp)(5, 12));
}
