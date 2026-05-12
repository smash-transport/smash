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
  T slope_;
  /// y-axis intercept of the linear interpolation.
  T yintercept_;
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
   *                           <tt>None</tt> and <tt>Constant</tt>.
   *
   * \return The interpolation function.
   * \throw std::out_of_range if values outside of the boundaries of the
   *                          underlying data are tried to be accessed and
   *                          extrapolation is disabled.
   * \throw std::invalid_argument if unsupported extrapolation type is
   *                              requested.
   */
  InterpolateDataLinear(
      const std::vector<T>& x, const std::vector<T>& y,
      const ExtrapolationType extrapolation_type = ExtrapolationType::None);

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
  std::vector<T> x_;
  /// Piecewise linear interpolation using f(x_i)
  std::vector<InterpolateLinear<T>> f_;
  /// Extrapolation type
  ExtrapolationType extrapolation_type_;
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

/**
 * Perform a trilinear 1st order interpolation
 *
 * Assume, we seek the value of a function f at position (x, y, z). We know the
 * position (x, y, z) lies within a 3D cube, for which the values of the
 * function f are known at each corner (f1, ..., f8). We can now interpolate
 * those values trilinearly to obtain an estimate of f at position (x, y, z).
 *
 * \param[in] ax fraction of the step in x-direction
 * \param[in] ay fraction of the step in y-direction
 * \param[in] az fraction of the step in z-direction
 * \param[in] f1 Value at the lower left front corner of the cube
 * \param[in] f2 Value at the lower right front corner of the cube
 * \param[in] f3 Value at the upper left front corner of the cube
 * \param[in] f4 Value at the upper right front corner of the cube
 * \param[in] f5 Value at the lower left back corner of the cube
 * \param[in] f6 Value at the lower right back corner of the cube
 * \param[in] f7 Value at the upper left back corner of the cube
 * \param[in] f8 Value at the upper right back corner of the cube
 *
 * - x-direction: lower left front to lower right front corner (f1 to f2)
 * - y-direction: lower left front to upper left front corner (f1 to f3)
 * - z-direction: lower left front to lower left back corner (f1 to f5)
 *
 * \return Interpolated value
 */
template <typename T>
T interpolate_trilinear(T ax, T ay, T az, T f1, T f2, T f3, T f4, T f5, T f6,
                        T f7, T f8) {
  assert(ax >= 0 && ax <= 1);
  assert(ay >= 0 && ay <= 1);
  assert(az >= 0 && az <= 1);
  T res = az * (ax * (ay * f8 + (1.0 - ay) * f6) +
                (1.0 - ax) * (ay * f7 + (1.0 - ay) * f5)) +
          (1 - az) * (ax * (ay * f4 + (1.0 - ay) * f2) +
                      (1.0 - ax) * (ay * f3 + (1.0 - ay) * f1));
  return res;
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
    const ExtrapolationType extrapolation_type) {
  extrapolation_type_ = extrapolation_type;
  assert(x.size() == y.size());
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
  double first_x = x_.front();
  double last_x = x_.back();
  if (x0 < first_x || x0 > last_x) {
    if (extrapolation_type_ == ExtrapolationType::None) {
      std::stringstream error_msg{};
      error_msg
          << "InterpolateDataLinear only accepts x values within the "
             "range of the underlying data when extrapolation is not specified."
          << "\nx value " << x0 << " is out of bounds.";
      throw std::out_of_range(error_msg.str());
    } else if (extrapolation_type_ == ExtrapolationType::Constant) {
      return (x0 < first_x) ? f_.front()(first_x) : f_.back()(last_x);
    } else {
      throw std::invalid_argument(
          "The provided extrapolation type is not supported. Valid types are "
          "'None' and 'Constant'.");
    }
  }
  // Find the piecewise linear interpolation corresponding to x0.
  size_t i = find_index(x_, x0);
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
   *                           <tt>None</tt> and <tt>Constant</tt>.
   *
   * \return The interpolation function.
   * \throw std::runtime_error if vectors x and y have different length.
   * \throw std::runtime_error if less than 3 data points are provided.
   * \throw std::out_of_range if values outside of the boundaries of the
   *                          underlying data are tried to be accessed and
   *                          extrapolation is disabled.
   * \throw std::invalid_argument if unsupported extrapolation type is
   *                              requested.
   */
  InterpolateDataSpline(
      const std::vector<double>& x, const std::vector<double>& y,
      const ExtrapolationType extrapolation_type = ExtrapolationType::None);

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
  ExtrapolationType extrapolation_type_;
  /// First x value of underlying data.
  double first_x_;
  /// Last x value of underlying data.
  double last_x_;
  /// First y value of underlying data.
  double first_y_;
  /// Last y value of underlying data.
  double last_y_;
  /// GSL iterator for interpolation lookups.
  gsl_interp_accel* acc_;
  /// GSL spline.
  gsl_spline* spline_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_INTERPOLATION_H_
