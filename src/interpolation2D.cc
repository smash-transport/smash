/*
 *    Copyright (c) 2020 -
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "smash/interpolation2D.h"

#include <initializer_list>
#include <iostream>

namespace smash {

InterpolateData2DSpline::InterpolateData2DSpline(const std::vector<double>& x,
                                                 const std::vector<double>& y,
                                                 const std::vector<double>& z) {
  const int M = x.size();
  const int N = y.size();

  if (z.size() != N * M) {
    throw std::runtime_error(
        "Dimensions not suitable for 2D interpolation. DIM(z) != DIM(x) * "
        "DIM(y).");
  }

  if (M < 4 || N < 4) {
    throw std::runtime_error(
        "Need at least 4 data points for bicubic spline interpolation.");
  }

  // Assign lower and upper bounds for constant extrapolation
  first_x_ = x.front();
  last_x_ = x.back();
  first_y_ = y.front();
  last_y_ = y.back();

  // cast vectors into arrays, as GSL functions can only handle arrays
  const double* xa = &x[0];
  const double* ya = &y[0];
  const double* za = &z[0];

  // Create accelerator objects (interpolation lookups)
  xacc_ = gsl_interp_accel_alloc();
  yacc_ = gsl_interp_accel_alloc();

  // Initialize bicubic spline interpolation
  spline_ = gsl_spline2d_alloc(gsl_interp2d_bicubic, M, N);
  gsl_spline2d_init(spline_, xa, ya, za, M, N);
}

InterpolateData2DSpline::~InterpolateData2DSpline() {
  gsl_spline2d_free(spline_);
  gsl_interp_accel_free(xacc_);
  gsl_interp_accel_free(yacc_);
}

double InterpolateData2DSpline::operator()(double xi, double yi) const {
  // constant extrapolation at the edges
  xi = (xi < first_x_) ? first_x_ : xi;
  xi = (xi > last_x_) ? last_x_ : xi;
  yi = (yi < first_y_) ? first_y_ : yi;
  yi = (yi > last_y_) ? last_y_ : yi;

  // bicubic spline interpolation
  return gsl_spline2d_eval(spline_, xi, yi, xacc_, yacc_);
}

}  // namespace smash
