/*
 *
 *    Copyright (c) 2012-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include "include/tabulationnd.h"

template class TabulationND<1>;
template class TabulationND<2>;

double TabulationND<1>::get_linear(double x) const {
    assert(x >= x0_ && x <= x1_);
  const double index_double = (x - x0_) * inv_dx_;
  // cast truncates, gets lower index
  const size_t lower = static_cast<size_t>(index_double);
  const size_t upper = std::min(lower + 1, n_ - 1);
  // TODO: Make more efficient. Looks stupid atm.
  const double x_lower = x0_ + lower * dx_;
  const double frac = (x - x_lower) * inv_dx_;
  return values_[lower] * (1 - frac) + values_[upper] * frac;
}

double TabulationND<2>::get_linear(const double x, const double y) const {
    assert(x >= x0_ && x <= x1_);
    assert(y >= y0_ && y <= y1_);
  // lower index

  const double x_idx_d = (x - x0_) * inv_dx_;
  const double y_idx_d = (y - y0_) * inv_dy_;

  const size_t x_idx = static_cast<size_t>(x_idx_d);
  const size_t y_idx = static_cast<size_t>(y_idx_d);

  // get fraction of x between x_lower, x_lower + dx etc.
  const double dx = x_idx_d - x_idx;
  const double dy = y_idx_d - y_idx;

  const double f_xy = values_[x_idx + y_idx * nx_] * (1 - dx) * (1 - dy) +
                      values_[x_idx + 1 + y_idx * nx_] * dx * (1 - dy) +
                      values_[x_idx + (y_idx + 1) * nx_] * (1 - dx) * dy +
                      values_[x_idx + 1 + (y_idx + 1) * nx_] * dx * dy;

  return f_xy;
}

double TabulationND<2>::get_closest(const double x, const double y) const {
  const double x_idx_d = (x - x0_) * inv_dx_;
  const double y_idx_d = (y - y0_) * inv_dy_;

  const size_t x_idx = static_cast<size_t>(x_idx_d);
  const size_t y_idx = static_cast<size_t>(y_idx_d);

  return values_[x_idx + nx_ * y_idx];
}
