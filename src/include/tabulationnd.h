/*
 *
 *    Copyright (c) 2012-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_TABULATIONND_H
#define SRC_INCLUDE_TABULATIONND_H

#include <cmath>
#include <functional>
#include <vector>

template <int dim>
class TabulationND {};

template <>
class TabulationND<1> {
  TabulationND(double x0, double x1, double dx, std::function<double(double)> f)
      : x0_{x0}, x1_{x1}, dx_{dx}, inv_dx_{1 / dx} {
    n_ = ceil((x1_ - x0_) * inv_dx_) + 1;
    values_.resize(n_);
    double x = x0;
    for (size_t i = 0; i < n_; i++, x += dx_) {
      values_[i] = f(x);
    }
  }

  double get_linear(double x);

  double step() const { return dx_; }
  double range() const { return x1_ - x0_; }
  double x0() const { return x0_; }
  double x1() const { return x1_; }
  size_t n_points() const { return n_; }

 private:
  const double x0_, x1_, dx_, inv_dx_;
  size_t n_;
  std::vector<double> values_;
};

double TabulationND<1>::get_linear(double x) {
  const double x_lower = (x - x0_) * inv_dx_;
  // cast truncates, gets lower index
  size_t index = static_cast<size_t>(x_lower);
  const double frac = (x - x_lower) * inv_dx_;
  return values_[index] * (1 - frac) + values_[index + 1] * frac;
}

#endif