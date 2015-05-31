/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "include/tabulation.h"

#include <gsl/gsl_integration.h>
#include <algorithm>

namespace Smash {

Tabulation::Tabulation(float x_min, float range, unsigned int N,
                       IntegParam ip, IntegrandFunction f)
                      : x_min_(x_min), dx_(range/N), ip_(ip), func_(f) {
  values_.resize(N);
  for (unsigned int i = 0; i < N; i++) {
    values_[i] = calculate_value(x_min_ + i*dx_);
  }
}

float Tabulation::calculate_value(float x) {
  gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(1000);
  ip_.srts = x;
  const gsl_function integrand = {func_, &ip_};
  const size_t subintervals_max = 500;
  const int gauss_points = 2;
  const double accuracy_absolute = 1.0e-5;
  const double accuracy_relative = 1.0e-3;
  double integral_value, integral_error;

  gsl_integration_qag(&integrand, ip_.type->minimum_mass(), ip_.srts - ip_.m2,
                      accuracy_absolute, accuracy_relative,
                      subintervals_max, gauss_points, workspace,
                      &integral_value, &integral_error);

  gsl_integration_workspace_free(workspace);

  return integral_value;
}

float Tabulation::get_value(float x) const {
  if (x < x_min_) {
    return 0.;
  }
  // look up tabulated values
  unsigned int n = static_cast<unsigned int>(std::round((x - x_min_) / dx_));
  if (n >= values_.size()) {
    return values_.back();
  } else {
    return values_[n];
  }
}

}  // namespace Smash
