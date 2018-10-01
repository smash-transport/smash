/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_INTERPOLATION_H_
#define SRC_INCLUDE_INTERPOLATION_H_

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace smash {

/**
 * Represent a linear interpolation.
 *
 * \param T Type of interpolated values.
 */
template <typename T>
class InterpolateLinear {
 public:
  /// Slope of the linear interpolation.
  T slope_;
  /// y-axis intercept of the linear interpolation.
  T yintercept_;
  /**
   * Linear interpolation given two points (x0, y0) and (x1, y1).
   *
   * \return The interpolation function.
   */
  InterpolateLinear(T x0, T y0, T x1, T y1);

  /**
   * Calculate spline interpolation at x.
   *
   *  \param x Interpolation argument.
   *  \return Interpolated value.
   */
  T operator()(T x) const;
};

/**
 * Represent a piecewise linear interpolation.
 *
 * \param T Type of interpolated values.
 */
template <typename T>
class InterpolateDataLinear {
 public:
  /**
   * Interpolate function f given discrete samples f(x_i) = y_i.
   *
   * \param x x-values.
   * \param y y-values.
   * \return The interpolation function.
   *
   * Piecewise linear interpolation is used.
   * Values outside the given samples will use the outmost linear
   * interpolation.
   */
  InterpolateDataLinear(const std::vector<T>& x, const std::vector<T>& y);
  /**
   * Calculate spline interpolation at x.
   *
   *  \param x Interpolation argument.
   *  \return Interpolated value.
   */
  T operator()(T x) const;

 private:
  /// x_i
  std::vector<T> x_;
  /// Piecewise linear interpolation using f(x_i)
  std::vector<InterpolateLinear<T>> f_;
};

template <typename T>
InterpolateLinear<T>::InterpolateLinear(T x0, T y0, T x1, T y1) {
  assert(x0 != x1);
  slope_ = (y1 - y0) / (x1 - x0);
  yintercept_ = y0 - slope_ * x0;
}

template <typename T>
T InterpolateLinear<T>::operator()(T x) const {
  return slope_ * x + yintercept_;
}

/// Represent a permutation.
using Permutation = std::vector<size_t>;

/**
 * Calculate the permutations necessary for sorting a vector.
 *
 * \tparam Cmp Type of comparison function.
 * \param v Vector to be sorted.
 * \param compare Comparison function (see `std::sort`).
 * \return Vector of indices into the original vector.
 */
template <typename T, typename Cmp>
Permutation generate_sort_permutation(std::vector<T> const& v, Cmp compare) {
  Permutation p(v.size());
  std::iota(p.begin(), p.end(), 0);
  std::sort(p.begin(), p.end(),
            [&](size_t i, size_t j) { return compare(v[i], v[j]); });
  return p;
}

/**
 * Apply a permutation to a vector.
 *
 * \tparam T Type of values to be permuted.
 * \param v Vector to be permuted.
 * \param p Permutation to be applied.
 * \return Permuted vector.
 */
template <typename T>
std::vector<T> apply_permutation(const std::vector<T>& v,
                                 const Permutation& p) {
  std::vector<T> copied_v = v;
  std::transform(p.begin(), p.end(), copied_v.begin(),
                 [&](size_t i) { return v[i]; });
  return copied_v;
}

/**
 * Check whether two components have the same value in a sorted vector x.
 *
 * Throws an exception if duplicates are encountered.
 *
 * \tparam T Type of values to be checked for duplicates.
 * \param x Vector to be checked for duplicates.
 * \param error_position String used in the error message, indicating where
 *                       the error originated.
 */
template <typename T>
void check_duplicates(const std::vector<T>& x,
                      const std::string& error_position) {
  auto it = std::adjacent_find(x.begin(), x.end());
  if (it != x.end()) {
    std::stringstream error_msg;
    error_msg << error_position << ": Each x value must be unique. \"" << *it
              << "\" was found twice.";
    throw std::runtime_error(error_msg.str());
  }
}

template <typename T>
InterpolateDataLinear<T>::InterpolateDataLinear(const std::vector<T>& x,
                                                const std::vector<T>& y) {
  assert(x.size() == y.size());
  const size_t n = x.size();
  const auto p = generate_sort_permutation(
      x, [&](T const& a, T const& b) { return a < b; });
  x_ = apply_permutation(x, p);
  check_duplicates(x_, "InterpolateDataLinear");
  std::vector<T> y_sorted = std::move(apply_permutation(y, p));
  f_.reserve(n - 1);
  for (size_t i = 0; i < n - 1; i++) {
    f_.emplace_back(
        InterpolateLinear<T>(x_[i], y_sorted[i], x_[i + 1], y_sorted[i + 1]));
  }
}

/**
 * Find the index in v that corresponds to the last value strictly smaller
 * than x. If no such value exists, the first value is returned.
 *
 * This assumes v is sorted and uses a binary search.
 *
 * \tparam T Type of values to be compared to x.
 * \param v Vector to be searched.
 * \param x Upper bound for indexed value.
 * \return Largest index corresponding to value below upper bound.
 *
 * Example:
 * >>> std::vector<int> x = { 0, 2, 4, 6, 8, 10 };
 * >>> find_index(x, 2)
 * 0
 * >>> find_index(x, 3)
 * 1
 */
template <typename T>
size_t find_index(const std::vector<T>& v, T x) {
  const auto it = std::lower_bound(v.begin(), v.end(), x);
  if (it == v.begin()) {
    return 0;
  } else {
    return it - 1 - v.begin();
  }
}

template <typename T>
T InterpolateDataLinear<T>::operator()(T x0) const {
  // Find the piecewise linear interpolation corresponding to x0.
  size_t i = find_index(x_, x0);
  if (i >= f_.size()) {
    // We don't have a linear interpolation beyond the last point in x_.
    // Use the last linear interpolation instead.
    i = f_.size() - 1;
  }
  return f_[i](x0);
}

/// Represent a cubic spline interpolation.
class InterpolateDataSpline {
 public:
  /**
   * Interpolate function f given discrete samples f(x_i) = y_i.
   *
   * \param x x-values.
   * \param y y-values.
   * \return The interpolation function.
   *
   * Cubic spline interpolation is used.
   * Values outside the given samples will use the outmost sample
   * as a constant extrapolation.
   */
  InterpolateDataSpline(const std::vector<double>& x,
                        const std::vector<double>& y);

  /// Destructor
  ~InterpolateDataSpline();

  /**
   * Calculate spline interpolation at x.
   *
   *  \param x Interpolation argument.
   *  \return Interpolated value.
   */
  double operator()(double x) const;

 private:
  /// First x value.
  double first_x_;
  /// Last x value.
  double last_x_;
  /// First y value.
  double first_y_;
  /// Last y value.
  double last_y_;
  /// GSL iterator for interpolation lookups.
  gsl_interp_accel* acc_;
  /// GSL spline.
  gsl_spline* spline_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_INTERPOLATION_H_
