/*
 *    Copyright (c) 2015-2018,2020,2026
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "smash/interpolation.h"

#include <iostream>

namespace smash {

InterpolateDataSpline::InterpolateDataSpline(
    const std::vector<double>& x, const std::vector<double>& y,
    const ExtrapolationType extrapolation_type) {
  extrapolation_type_ = extrapolation_type;
  switch (extrapolation_type_) {
    case ExtrapolationType::None:
    case ExtrapolationType::Zero:
    case ExtrapolationType::Constant:
      break;
    default:
      throw std::invalid_argument(
          "The provided extrapolation type is not supported. Valid types are "
          "'None', 'Zero', and 'Constant'.");
  }
  const auto N = x.size();
  if (y.size() != N) {
    throw std::invalid_argument(
        "The interpolation requires two vectors of equal length.");
  }
  if (N < 3) {
    throw std::invalid_argument(
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
  if (xi < first_x_ || xi > last_x_) {
    if (extrapolation_type_ == ExtrapolationType::None) {
      std::stringstream error_msg{
          "InterpolateDataSpline only accepts x values within the range of the "
          "underlying data when an extrapolation type is not specified."};
      error_msg << "\nx value " << xi << " is out of bounds.";
      throw std::out_of_range(error_msg.str());
    } else if (extrapolation_type_ == ExtrapolationType::Zero) {
      return 0.;
    } else if (extrapolation_type_ == ExtrapolationType::Constant) {
      return (xi < first_x_) ? first_y_ : last_y_;
    }
  }
  // cubic spline interpolation
  return gsl_spline_eval(spline_, xi, acc_);
}

}  // namespace smash
