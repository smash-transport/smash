/*
 *
 *    Copyright (c) 2015-2018,2020,2022,2024,2026
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_INTERPOLATION_H_
#define SRC_INCLUDE_SMASH_INTERPOLATION_H_

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"

#include "smash/constants.h"
#include "smash/forwarddeclarations.h"

namespace smash {

/**
 * Represent a linear interpolation.
 *
 * \param T Type of interpolated values.
 */
template <typename T>
class InterpolateLinear {
 public:
  /**
   * Linear interpolation given two points (x0, y0) and (x1, y1).
   *
   * \return The interpolation function.
   */
  InterpolateLinear(T x0, T y0, T x1, T y1);

  /**
   * Calculate linear interpolation at x.
   *
   * \param x Interpolation argument.
   * \return Interpolated value.
   */
  T operator()(T x) const;

 private:
  /// Slope of the linear interpolation.
  T slope_{};
  /// y-axis intercept of the linear interpolation.
  T yintercept_{};
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
   * Piecewise linear interpolation is used.
   *
   * \param x x-values.
   * \param y y-values.
   * \param extrapolation_type Type of extrapolation for requested x_i values
   *                           that are out of bounds. Extrapolation is by
   *                           default disabled. Possible types are
   *                           <tt>None</tt>, <tt>Zero</tt>, <tt>Constant</tt>,
   *                           and <tt>Linear</tt>.
   *
   * \return The interpolation function.
   * \throw std::invalid_argument if vectors x and y have different length.
   * \throw std::out_of_range if values outside of the boundaries of the
   *                          underlying data are tried to be accessed and
   *                          extrapolation is disabled.
   * \throw std::invalid_argument if unsupported extrapolation type is
   *                              requested.
   */
  InterpolateDataLinear(
      const std::vector<T>& x, const std::vector<T>& y,
      ExtrapolationType extrapolation_type = ExtrapolationType::None);

  /**
   * Calculate linear interpolation at x.
   *
   * \param x Interpolation argument.
   * \return Interpolated value.
   *
   * \throw std::out_of_range if values outside of the boundaries of the
   *                          underlying data are tried to be accessed and
   *                          extrapolation is disabled.
   */
  T operator()(T x) const;

 private:
  /// x_i
  std::vector<T> x_{};
  /// Piecewise linear interpolation using f(x_i)
  std::vector<InterpolateLinear<T>> f_{};
  /// Extrapolation type
  ExtrapolationType extrapolation_type_ = ExtrapolationType::None;
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
 * \tparam T Type of values to be checked for duplicates.
 * \param x Vector to be checked for duplicates.
 * \param error_position String used in the error message, indicating where the
 *                       error originated.
 *
 * \throw std::runtime_error if duplicates are encountered.
 */
template <typename T>
void check_duplicates(const std::vector<T>& x,
                      const std::string& error_position) {
  auto it = std::adjacent_find(x.begin(), x.end());
  if (it != x.end()) {
    std::stringstream error_msg{};
    error_msg << error_position << ": Each x value must be unique. \"" << *it
              << "\" was found twice.";
    throw std::runtime_error(error_msg.str());
  }
}

template <typename T>
InterpolateDataLinear<T>::InterpolateDataLinear(
    const std::vector<T>& x, const std::vector<T>& y,
    const ExtrapolationType extrapolation_type)
    : extrapolation_type_{extrapolation_type} {
  switch (extrapolation_type_) {
    case ExtrapolationType::None:
    case ExtrapolationType::Zero:
    case ExtrapolationType::Constant:
    case ExtrapolationType::Linear:
      break;
    default:
      throw std::invalid_argument(
          "The provided extrapolation type is not supported. Valid types are "
          "'None', 'Zero', 'Constant', and 'Linear'.");
  }
  if (x.size() != y.size()) {
    throw std::invalid_argument(
        "The interpolation requires two vectors of equal length.");
  }
  const size_t n = x.size();
  const auto p = generate_sort_permutation(
      x, [&](T const& a, T const& b) { return a < b; });
  x_ = apply_permutation(x, p);
  check_duplicates(x_, "InterpolateDataLinear");
  std::vector<T> y_sorted = apply_permutation(y, p);
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
  const double first_x = x_.front();
  const double last_x = x_.back();
  if (x0 < first_x || x0 > last_x) {
    if (extrapolation_type_ == ExtrapolationType::None) {
      std::ostringstream error_msg{
          "InterpolateDataLinear only accepts x values within the range of the "
          "underlying data\nwhen an extrapolation type is not specified. ",
          std::ios::ate};
      error_msg << "x value " << x0 << " is out of bounds.";
      throw std::out_of_range(error_msg.str());
    } else if (extrapolation_type_ == ExtrapolationType::Zero) {
      return 0.;
    } else if (extrapolation_type_ == ExtrapolationType::Constant) {
      return (x0 < first_x) ? f_.front()(first_x) : f_.back()(last_x);
    } else if (extrapolation_type_ == ExtrapolationType::Linear) {
      return (x0 < first_x) ? f_.front()(x0) : f_.back()(x0);
    }
  }
  // Find the piecewise linear interpolation corresponding to x0.
  const size_t i = find_index(x_, x0);
  return f_[i](x0);
}

/// Represent a cubic spline interpolation.
class InterpolateDataSpline {
 public:
  /**
   * Interpolate function f given discrete samples f(x_i) = y_i.
   * Cubic spline interpolation is used.
   *
   * \param x x-values.
   * \param y y-values.
   * \param extrapolation_type Type of extrapolation for requested x_i values
   *                           that are out of bounds. Extrapolation is by
   *                           default disabled. Possible types are
   *                           <tt>None</tt>, <tt>Zero</tt>, and
   *                           <tt>Constant</tt>.
   *
   * \return The interpolation function.
   * \throw std::invalid_argument if vectors x and y have different length.
   * \throw std::invalid_argument if less than 3 data points are provided.
   * \throw std::out_of_range if values outside of the boundaries of the
   *                          underlying data are tried to be accessed and
   *                          extrapolation is disabled.
   * \throw std::invalid_argument if unsupported extrapolation type is
   *                              requested.
   */
  InterpolateDataSpline(
      const std::vector<double>& x, const std::vector<double>& y,
      ExtrapolationType extrapolation_type = ExtrapolationType::None);

  /// Destructor
  ~InterpolateDataSpline();

  /**
   * Calculate spline interpolation at x.
   *
   * \param x Interpolation argument.
   * \return Interpolated value.
   *
   * \throw std::out_of_range if values outside of the boundaries of the
   *                          underlying data are tried to be accessed and
   *                          extrapolation is disabled.
   */
  double operator()(double x) const;

 private:
  /// Extrapolation type.
  ExtrapolationType extrapolation_type_ = ExtrapolationType::None;
  /// First x value of underlying data.
  double first_x_ = smash_NaN<double>;
  /// Last x value of underlying data.
  double last_x_ = smash_NaN<double>;
  /// First y value of underlying data.
  double first_y_ = smash_NaN<double>;
  /// Last y value of underlying data.
  double last_y_ = smash_NaN<double>;
  /// GSL iterator for interpolation lookups.
  gsl_interp_accel* acc_ = nullptr;
  /// GSL spline.
  gsl_spline* spline_ = nullptr;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_INTERPOLATION_H_
