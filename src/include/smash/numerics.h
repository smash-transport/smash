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
#include <cassert>
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

namespace smash {

namespace detail {

/**
 * Compare whether two floating-point numbers are approximately equal à la Knuth
 * up to a given tolerance. On top of Knuth's tolerance predicate some corner
 * cases are treated and the caller can specify a threshold as last parameter to
 * make the test consider numbers equal if the absolute value of their
 * difference is below of it.
 *
 * \param[in] x First of the two numbers.
 * \param[in] y Second of the two numbers.
 * \param[in] epsilon The relative tolerance for the test.
 * \param[in] threshold Threshold for the number comparison. By default this is
 * zero, implying no threshold is considered.
 *
 * \return \c false if either \c x or \c y is not a finite number, provided that
 * the type supports non-numeric representations (i.e. is infinite or NAN);
 * \return \c true if <tt>x == y</tt>;
 * \return \c true if <tt>x == 0</tt> and if \f$ |x| \le \varepsilon\f$;
 * \return \c true if <tt>y == 0</tt> and if \f$ |y| \le \varepsilon\f$;
 * \return \c true if \f$ |x - y| \le M_\mathrm{threshold} \f$;
 * \return \c true if \f$ |x - y| \le \varepsilon \cdot \max(|x|, |y|) \f$
 * (Knuth's tolerance predicate);
 * \return \c false otherwise.
 */
template <typename N, typename = std::enable_if_t<std::is_floating_point_v<N>>>
bool almost_equal_knuthish(const N x, const N y, const N epsilon,
                           const N threshold = N{0.0}) noexcept {
  assert(epsilon > 0);
  assert(threshold >= 0);
  if constexpr (std::numeric_limits<N>::is_iec559) {
    if (!std::isfinite(x) || !std::isfinite(y)) {
      return false;
    }
  }
  if (x == y)
    return true;
  else if (x == 0)
    return std::abs(y) <= epsilon;
  else if (y == 0)
    return std::abs(x) <= epsilon;
  else
    return std::abs(x - y) <= threshold ||
           std::abs(x - y) <= epsilon * std::max(std::abs(x), std::abs(y));
}

}  // namespace detail

/**
 * Checks whether two floating-point numbers are almost equal. This is done
 * using \c smash::really_small as relative tolerance. All numbers are tested
 * for equality, no matter which order of magnitude they have.
 *
 * \see detail::almost_equal_knuthish for more information.
 */
template <typename N, typename = std::enable_if_t<std::is_floating_point_v<N>>>
bool almost_equal(const N x, const N y) {
  return detail::almost_equal_knuthish<N>(x, y, static_cast<N>(really_small));
}

/**
 * Like \c smash::almost_equal, but using a less strict tolerance,
 * <tt>smash::small_number</tt>. Furthermore, numbers smaller than a given
 * threshold, <tt>smash::really_small</tt>, are now considered equal (in the
 * sense that, for SMASH physics, their difference has no physical effect).
 */
template <typename N, typename = std::enable_if_t<std::is_floating_point_v<N>>>
bool almost_equal_physics(const N x, const N y) {
  const auto threshold = static_cast<N>(really_small);
  const auto epsilon = static_cast<N>(small_number);
  return detail::almost_equal_knuthish<N>(x, y, epsilon, threshold);
}

/**
 * \brief Returns whether any element in a collection is NaN.
 *
 * This function iterates through the elements of a collection and checks if any
 * of them is NaN using \c std::isnan. NaN is a special floating-point value
 * that represents undefined or unrepresentable values.
 *
 * \tparam T Iterable container of numeric values (defaults to \c
 * std::initializer_list<double>).
 * \param[in] collection The collection to be checked for NaN values.
 *
 * \return \c true if any element in the collection is NaN,
 * \return \c false otherwise.
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
