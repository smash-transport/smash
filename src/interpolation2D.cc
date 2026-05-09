/*
 *    Copyright (c) 2020,2022,2026
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "smash/interpolation2D.h"

#include <algorithm>
#include <initializer_list>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace smash {

InterpolateData2DSpline::InterpolateData2DSpline(
    const std::vector<double>& x, const std::vector<double>& y,
    const std::vector<double>& z, const ExtrapolationType extrapolation_type) {
  extrapolation_type_ = extrapolation_type;
  const size_t M = x.size();
  const size_t N = y.size();

  if (z.size() != N * M) {
    throw std::runtime_error(
        "Dimensions not suitable for 2D interpolation. DIM(z) != DIM(x) * "
        "DIM(y).");
  }

  if (M < 4 || N < 4) {
    throw std::runtime_error(
        "Need at least 4 data points in each dimension for bicubic spline "
        "interpolation.");
  }

  if (!std::is_sorted(x.begin(), x.end()) ||
      !std::is_sorted(y.begin(), y.end())) {
    throw std::runtime_error(
        "x and y values must be strictly increasing, i.e. the vectors have to "
        "be sorted by size of their values. This is required by GSL.");
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
  if ((xi < first_x_ || xi > last_x_) || (yi < first_y_ || yi > last_y_)) {
    if (extrapolation_type_ == ExtrapolationType::None) {
      std::stringstream error_msg{};
      error_msg << "InterpolateData2DSpline only accepts x and y values within "
                   "the range of the underlying data when extrapolation is not "
                   "specified."
                << "\nx value " << xi << " or y value " << yi
                << " are out of bounds.";
      throw std::out_of_range(error_msg.str());
    } else if (extrapolation_type_ == ExtrapolationType::Constant_value) {
      // constant extrapolation at the edges
      xi = (xi < first_x_) ? first_x_ : xi;
      xi = (xi > last_x_) ? last_x_ : xi;
      yi = (yi < first_y_) ? first_y_ : yi;
      yi = (yi > last_y_) ? last_y_ : yi;
    }
  }

  // bicubic spline interpolation
  return gsl_spline2d_eval(spline_, xi, yi, xacc_, yacc_);
}

}  // namespace smash
