/*
 *
 *    Copyright (c) 2014-2018,2020,2022-2023
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_NUMERICS_H_
#define SRC_INCLUDE_SMASH_NUMERICS_H_

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <type_traits>

#include "constants.h"

/**
 * \file
 *
 * Generic numerical functions
 *
 * This file collects generic numerical functions such as for
 * approximate equality checks between two floating point values.
 */

namespace details {

/**
 * \brief Relative (by default) approximate equality for floating-point numbers.
 *
 * Uses the **maximum magnitude** as scale:
 * \f[
 *   |a - b| \le \text{rel\_tol} \cdot \max(|a|, |b|)
 * \f]
 *
 * If \p abs_tol > 0, uses a mixed criterion:
 * \f[
 *   |a - b| \le \max\!\Big(\text{abs\_tol},\ \text{rel\_tol} \cdot \max(|a|,
 * |b|)\Big). \f]
 *
 * - Fast-path equality returns true (also handles +0 == -0 and equal
 * infinities).
 * - If either operand is non-finite (NaN or infinity) and not equal by the fast
 * path, returns false.
 */
template <typename Float,
          typename = std::enable_if_t<std::is_floating_point_v<Float>>>
bool almost_equal_knuthish(Float a, Float b, Float rel_tol,
                           Float abs_tol = Float(0)) noexcept {
  if (a == b)
    return true;

  if constexpr (std::numeric_limits<Float>::is_iec559) {
    if (!std::isfinite(a) || !std::isfinite(b))
      return false;
  }

  const Float a_abs = std::abs(a);
  const Float b_abs = std::abs(b);
  const Float diff = std::abs(a - b);
  const Float scale = std::max(a_abs, b_abs);
  const Float rlimit = rel_tol * scale;

  if (abs_tol > Float(0)) {
    return diff <= std::max(abs_tol, rlimit);
  }
  return diff <= rlimit;
}

}  // namespace details

namespace smash {

/**
 * \brief Checks two numbers for **relative-only** approximate equality.
 *
 * Pure relative criterion using \c smash::really_small and the max-magnitude
 * scale: \f[ |x - y| \le \texttt{smash::really\_small} \cdot \max(|x|, |y|).
 * \f]
 *
 * \note No absolute floor is used.
 */
template <typename Float,
          typename = std::enable_if_t<std::is_floating_point_v<Float>>>
bool almost_equal(Float x, Float y) {
  return details::almost_equal_knuthish<Float>(
      x, y, static_cast<Float>(really_small), static_cast<Float>(0));
}

/**
 * \brief Like smash::almost_equal but with an **absolute floor** = \c
 * smash::small_number.
 *
 * Mixed criterion with max-magnitude scaling:
 * \f[
 *   |x - y| \le \max\!\Big(\texttt{smash::small\_number},\
 *                           \texttt{smash::small\_number} \cdot \max(|x|,
 * |y|)\Big). \f]
 */
template <typename Float,
          typename = std::enable_if_t<std::is_floating_point_v<Float>>>
bool almost_equal_physics(Float x, Float y) {
  const auto eps = static_cast<Float>(small_number);
  return details::almost_equal_knuthish<Float>(x, y, eps, eps);
}

/**
 * \brief Returns whether any element in a collection is NaN.
 *
 * This function iterates through the elements of a collection and checks if any
 * of them is NaN using \c std::isnan. NaN is a special floating-point value
 * that represents undefined or unrepresentable values.
 *
 * \tparam T Iterable container of numeric values (defaults to \c
 * std::initializer_list<double>). \param collection The collection to be
 * checked for NaN values. \return \c true if any element in the collection is
 * NaN, \c false otherwise.
 */
template <typename T = std::initializer_list<double>>
bool is_any_nan(const T& collection) {
  for (const auto& number : collection) {
    if (unlikely(std::isnan(number)))
      return true;
  }
  return false;
}

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_NUMERICS_H_
