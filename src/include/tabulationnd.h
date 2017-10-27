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
#include <iostream>
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
  double get_from_index(size_t i) const { return values_.at(i); }
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



template <>
class TabulationND<2> {
 public:
  TabulationND(const double x0, const double x1, const double y0,
               const double y1, const double dx, const double dy,
               std::function<double(double, double)> f)
      : x0_{x0},
        x1_{x1},
        y0_{y0},
        y1_{y1},
        dx_{dx},
        dy_{dy},
        inv_dx_{1 / dx},
        inv_dy_{1 / dy} {
    nx_ = std::ceil((x1_ - x0_) * inv_dx_ + 1);
    ny_ = std::ceil((y1_ - y0_) * inv_dy_ + 1);
    n_ = nx_ * ny_;

    std::cout << "nx, ny " << nx_ << " " << ny_ << "\n";

    values_.resize(n_);
    double x, y;
    for (size_t i = 0; i < ny_; i++) {
      for (size_t j = 0; j < nx_; j++) {
        x = x0_ + j * dx_;
        y = y0_ + i * dy_;
        values_[j + i * nx_] = f(x, y);
      }
    }
  }

  double get_linear(const double x, const double y) const;
  double get_closest(const double x, const double y) const;

 private:
  const double x0_, x1_, y0_, y1_, dx_, dy_, inv_dx_, inv_dy_;
  size_t n_, nx_, ny_;
  std::vector<double> values_;
};


#endif