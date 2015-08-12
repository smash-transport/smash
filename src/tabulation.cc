/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "include/tabulation.h"

namespace Smash {

Tabulation::Tabulation(float x_min, float range, int num,
                       std::function<double(float)> f)
                      : x_min_(x_min), inv_dx_(num/range) {
  values_.resize(num+1);
  const float dx = range/num;
  for (int i = 0; i <= num; i++) {
    values_[i] = f(x_min_ + i*dx);
  }
}

float Tabulation::get_value_step(float x) const {
  if (x < x_min_) {
    return 0.;
  }
  // this rounds correctly because float -> int conversion truncates
  const unsigned int n = (x - x_min_) * inv_dx_ + 0.5f;
  if (n >= values_.size()) {
    return values_.back();
  } else {
    return values_[n];
  }
}

float Tabulation::get_value_linear(float x) const {
  if (x < x_min_) {
    return 0.;
  }
  // here n is the lower index
  const unsigned int n = std::min((x - x_min_) * inv_dx_, values_.size() - 2.f);
  const float r = (x - x_min_) * inv_dx_ - n;
  return values_[n] * (1. - r) + values_[n + 1] * r;
}

}  // namespace Smash
