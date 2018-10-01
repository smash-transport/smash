/*
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "smash/tabulation.h"

namespace smash {

Tabulation::Tabulation(double x_min, double range, int num,
                       std::function<double(double)> f)
    : x_min_(x_min), x_max_(x_min + range), inv_dx_(num / range) {
  values_.resize(num + 1);
  const double dx = range / num;
  for (int i = 0; i <= num; i++) {
    values_[i] = f(x_min_ + i * dx);
  }
}

double Tabulation::get_value_step(double x) const {
  if (x < x_min_) {
    return 0.;
  }
  // this rounds correctly because double -> int conversion truncates
  const unsigned int n = (x - x_min_) * inv_dx_ + 0.5;
  if (n >= values_.size()) {
    return values_.back();
  } else {
    return values_[n];
  }
}

double Tabulation::get_value_linear(double x, Extrapolation extrapol) const {
  if (x < x_min_) {
    return 0.;
  }
  if (extrapol == Extrapolation::Zero && x > x_max_) {
    return 0.0;
  }
  if (extrapol == Extrapolation::Const && x > x_max_) {
    return values_.back();
  }
  const double index_double = (x - x_min_) * inv_dx_;
  // here n is the lower index
  const size_t n =
      std::min(static_cast<size_t>(index_double), values_.size() - 2);
  const double r = index_double - n;
  return values_[n] + (values_[n + 1] - values_[n]) * r;
}

}  // namespace smash
