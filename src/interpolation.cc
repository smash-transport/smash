/*
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "smash/interpolation.h"

#include <iostream>

namespace smash {

InterpolateDataSpline::InterpolateDataSpline(const std::vector<double>& x,
                                             const std::vector<double>& y) {
  const auto N = x.size();
  if (y.size() != N) {
    throw std::runtime_error(
        "Need two vectors of equal length for interpolation.");
  }
  if (N < 3) {
    throw std::runtime_error(
        "Need at least 3 data points for cubic spline interpolation.");
  }
  const auto p = generate_sort_permutation(
      x, [&](double const& a, double const& b) { return a < b; });
  const std::vector<double> sorted_x = apply_permutation(x, p);
  const std::vector<double> sorted_y = apply_permutation(y, p);
  check_duplicates(sorted_x, "InterpolateDataSpline");

  first_x_ = sorted_x.front();
  last_x_ = sorted_x.back();
  first_y_ = sorted_y.front();
  last_y_ = sorted_y.back();
  acc_ = gsl_interp_accel_alloc();
  spline_ = gsl_spline_alloc(gsl_interp_cspline, N);
  gsl_spline_init(spline_, &(*sorted_x.begin()), &(*sorted_y.begin()), N);
}

InterpolateDataSpline::~InterpolateDataSpline() {
  gsl_spline_free(spline_);
  gsl_interp_accel_free(acc_);
}

double InterpolateDataSpline::operator()(double xi) const {
  // constant extrapolation
  if (xi < first_x_) {
    return first_y_;
  }
  if (xi > last_x_) {
    return last_y_;
  }
  // cubic spline interpolation
  return gsl_spline_eval(spline_, xi, acc_);
}

}  // namespace smash
