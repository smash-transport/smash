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
 public:
  TabulationND(double x0, double x1, double dx, std::function<double(double)> f)
      : x0_{x0}, x1_{x1}, dx_{dx}, inv_dx_{1 / dx} {
    n_ = ceil((x1_ - x0_) * inv_dx_) + 1;
    values_.resize(n_);
    double x = x0;
    for (size_t i = 0; i < n_; i++, x += dx_) {
      values_[i] = f(x);
    }
  }

  // interpolate linear between tabulated values
  double get_linear(double x) const;
  double get_from_index(size_t i) const { return values_.at(i);}
  // getter functions
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

double TabulationND<1>::get_linear(double x) const {
  const double index_double = (x - x0_) * inv_dx_;
  // cast truncates, gets lower index
  const size_t lower = static_cast<size_t>(index_double);
  const size_t upper = std::min(lower + 1, n_- 1);
  // TODO: Make more efficient. Looks stupid atm.
  const double x_lower = x0_ + lower * dx_;
  const double frac = (x - x_lower) * inv_dx_;
  return values_[lower] * (1 - frac) + values_[upper] * frac;
}

#endif