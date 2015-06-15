/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_INTERPOLATION_H_
#define SRC_INCLUDE_INTERPOLATION_H_

#include <vector>
#include <cassert>
#include <cstddef>
#include <algorithm>

template <typename T>
class InterpolateLinear {
 public:
  T slope_, yintercept_;
  /// Linear interpolation given two points (x0, y0) and (x1, y1).
  ///
  /// Returns the interpolation function.
  InterpolateLinear(T x0, T y0, T x1, T y1);
  /// Calculate linear interpolation at x.
  T operator()(T x) const;
};

template <typename T>
class InterpolateData {
 public:
  /// Interpolate function f given discrete samples f(x_i) = y_i.
  //
  /// Returns the interpolation function.
  ///
  /// Piecewise linear interpolation is used.
  /// Values outside the given samples will use the outmost linear
  /// interpolation.
  InterpolateData(const std::vector<T>& x, const std::vector<T>& y);
  /// Calculate interpolation of f at x.
  T operator()(T x) const;

 private:
  std::vector<T> x_;
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

using Permutation = std::vector<size_t>;

/// Calculate the permutations necessary for sorting a vector.
template <typename T, typename Cmp>
Permutation generate_sort_permutation(std::vector<T> const& v, Cmp compare) {
  Permutation p(v.size());
  std::iota(p.begin(), p.end(), 0);
  std::sort(p.begin(), p.end(),
            [&](size_t i, size_t j) { return compare(v[i], v[j]); });
  return p;
}

/// Apply a permutation to a vector in place.
template <typename T>
void apply_permutation(std::vector<T>& v, const Permutation& p) {
  std::transform(p.begin(), p.end(), v.begin(), [&](size_t i) { return v[i]; });
}

template <typename T>
InterpolateData<T>::InterpolateData(const std::vector<T>& x,
                                    const std::vector<T>& y) {
  assert(x.size() == y.size());
  const size_t n = x.size();
  x_ = x;
  std::vector<T> y_sorted = y;
  const auto p = generate_sort_permutation(
      x, [&](T const& a, T const& b) { return a < b; });
  apply_permutation(x_, p);
  apply_permutation(y_sorted, p);
  f_.reserve(n - 1);
  for (size_t i = 0; i < n - 1; i++) {
    f_.emplace_back(
        InterpolateLinear<T>(x[i], y_sorted[i], x[i + 1], y_sorted[i + 1]));
  }
}

template <typename T>
T InterpolateData<T>::operator()(T x0) const {
  // Find the piecewise linear interpolation corresponding to x0.
  size_t i = 0;
  while (x_[i] <= x0 && i < f_.size()) {
    i++;
  }
  if (i != 0) {
    i--;
  }
  return f_[i](x0);
}

#endif  // SRC_INCLUDE_INTERPOLATION_H_
