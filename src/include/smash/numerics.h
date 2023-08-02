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

#include <cmath>
#include <initializer_list>

#include "constants.h"

/**
 * \file
 *
 * Generic numerical functions
 *
 * This file collects generic numerical functions such as for
 * approximate equality checks between to floating point values.
 */

namespace smash {
/**
 * Checks two numbers for relative approximate equality.
 *
 * \tparam N Number type.
 * \param x Left-hand side.
 * \param y Right-hand side.
 * \return true if the difference between x and y is less than or equal to
 * \f$\delta = \f$smash::really_small or that times the average of
 * \f$|x|\f$ and \f$|y|\f$:
 *
 * \f[|x - y| \stackrel{?}{\le} \delta \mbox{ or } |x - y|
 * \stackrel{?}{\le} \frac{|x| + |y|}{2} \cdot \delta\f]
 *
 * \see smash::really_small
 */
template <typename N>
bool almost_equal(const N x, const N y) {
  return (std::abs(x - y) <= N(really_small) ||
          std::abs(x - y) <=
              N(0.5 * really_small) * (std::abs(x) + std::abs(y)));
}
/**
 * Same as smash::almost_equal, but for physical checks like energy-momentum
 * conservation small_number is enough precision-wise
 *
 * \tparam N Number type.
 * \param x Left-hand side.
 * \param y Right-hand side.
 * \return true if the difference between x and y is less than or equal to
 * \f$\delta = \f$smash::small_number or that times the average of
 * \f$|x|\f$ and \f$|y|\f$:
 *
 * \see smash::small_number
 * \see smash::almost_equal
 */
template <typename N>
bool almost_equal_physics(const N x, const N y) {
  return (std::abs(x - y) <= N(small_number) ||
          std::abs(x - y) <=
              N(0.5 * small_number) * (std::abs(x) + std::abs(y)));
}

/**
 * This function iterates through the elements of a collection and checks if any
 * of them is NaN using the \c std::isnan function. NaN is a special
 * floating-point value that represents undefined or unrepresentable values.
 *
 * \tparam T The type of the collection. It can be any iterable container of
 *         numeric values.
 * \param collection The collection to be checked for NaN
 *        values.
 * \return \c true if any element in the collection is NaN, \c false otherwise
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
