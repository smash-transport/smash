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

static void swrite(std::ofstream& stream, double x) {
  stream.write(reinterpret_cast<const char*>(&x), sizeof(x));
}

static void swrite(std::ofstream& stream, size_t x) {
  // We want to support 32-bit and 64-bit platforms, so we store a 64-bit
  // integer on all platforms.
  const auto const_size_x = static_cast<uint64_t>(x);
  stream.write(reinterpret_cast<const char*>(&const_size_x), sizeof(const_size_x));
}

static void swrite(std::ofstream& stream, const std::vector<double> x) {
  swrite(stream, x.size());
  if (x.size() > 0) {
    stream.write(reinterpret_cast<const char*>(x.data()), sizeof(x[0]) * x.size());
  }
}

static void swrite(std::ofstream& stream, sha256::Hash x) {
  // The size is always the same, so there is no need to write it.
  stream.write(reinterpret_cast<const char*>(x.data()), sizeof(x[0]) * x.size());
}

void Tabulation::write(std::ofstream& stream, sha256::Hash hash) {
  swrite(stream, hash);
  swrite(stream, x_min_);
  swrite(stream, x_max_);
  swrite(stream, inv_dx_);
  swrite(stream, values_);
}

}  // namespace smash
