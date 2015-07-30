#include "include/interpolation.h"

#include <iostream>

InterpolateDataSpline::InterpolateDataSpline(const std::vector<double>& x,
                                             const std::vector<double>& y)
    : first_x_(x.front()), last_x_(x.back()), first_y_(y.front()), last_y_(y.back()) {
    const auto N = x.size();
    if (y.size() != N) {
      throw std::runtime_error("Need two vectors of equal length for interpolation.");
    }
    if (N < 3) {
      throw std::runtime_error("Need at least 3 data points for cubic spline interpolation.");
    }
    acc_ = gsl_interp_accel_alloc();
    spline_ = gsl_spline_alloc(gsl_interp_cspline, N);
    gsl_spline_init(spline_, &(*x.begin()), &(*y.begin()), N);
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
