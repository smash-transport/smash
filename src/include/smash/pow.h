/*
 *
 *    Copyright (c) 2016-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_POW_H_
#define SRC_INCLUDE_POW_H_

namespace smash {

/**
 * Efficient template for calculating integer powers using squaring.
 * \tparam T Type that implements multiplication.
 * \param[in] base
 * \param[in] exponent
 * \return base^exponent
 */
template <class T>
inline constexpr T pow_int(const T base, unsigned const exponent) {
  return (exponent == 0)
             ? 1
             : (exponent % 2 == 0)
                   ? pow_int(base, exponent / 2) * pow_int(base, exponent / 2)
                   : base * pow_int(base, exponent - 1);
}

/**
 * Efficient template for calculating the square.
 * \tparam T Type that implements multiplication.
 * \param[in] base value that gets squared
 * \return result of squaring base
 */
template <class T>
inline constexpr T square(const T base) {
  return pow_int(base, 2);
}

}  // namespace smash

#endif  // SRC_INCLUDE_POW_H_
