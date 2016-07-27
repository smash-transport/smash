/*
 *
 *    Copyright (c) 2016
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

/// Efficient templates for calculating integer powers.
#ifndef SRC_INCLUDE_POW_H_
#define SRC_INCLUDE_POW_H_

namespace Smash {

/// Calculate integer powers using squaring.
template<class T>
inline constexpr T pow(const T base, unsigned const exponent) {
    return (exponent == 0)     ? 1 :
           (exponent % 2 == 0) ? pow(base, exponent/2) * pow(base, exponent/2) :
                                 base * pow(base, exponent - 1);
}

template<class T>
inline constexpr T square(const T base) {
    return pow(base, 2);
}

}  // namespace Smash

#endif  // SRC_INCLUDE_POW_H_
