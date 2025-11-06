/*
 *
 *    Copyright (c) 2016-2018,2020,2022,2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_POW_H_
#define SRC_INCLUDE_SMASH_POW_H_

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
  return (exponent == 0) ? 1
         : (exponent % 2 == 0)
             ? pow_int(base, exponent / 2) * pow_int(base, exponent / 2)
             : base * pow_int(base, exponent - 1);
}

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_POW_H_
