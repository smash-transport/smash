/*
 *
 *    Copyright (c) 2014-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_NUMERICS_H_
#define SRC_INCLUDE_NUMERICS_H_

#include <cmath>
#include "constants.h"

/**
 * \file
 *
 * Generic numerical functions
 *
 * This file collects generic numerical functions such as for
 * approximate equality checks between to floating point values.
 *
 */

namespace smash {
/** checks two numbers for relative approximate equality.
 *
 * \param x
 * \param y
 *
 * \returns true the difference between x and y is less than or equal to
 * \f$\delta = \f$SMASH::really_small or that times the average of
 * \f$|x|\f$ and \f$|y|\f$:
 *
 * \f[|x - y| \stackrel{?}{\le} \delta \mbox{ or } |x - y|
 * \stackrel{?}{\le} \frac{|x| + |y|}{2} \cdot \delta\f]
 */
template <typename N>
bool almost_equal(const N x, const N y) {
  return (std::abs(x - y) <= N(really_small) ||
          std::abs(x - y) <=
              N(0.5 * really_small) * (std::abs(x) + std::abs(y)));
}
/** Same as above, but for physical checks like energy-momentum conser-
 * vation small_number is enough precision-wise
 */
template <typename N>
bool almost_equal_physics(const N x, const N y) {
  return (std::abs(x - y) <= N(small_number) ||
          std::abs(x - y) <=
              N(0.5 * small_number) * (std::abs(x) + std::abs(y)));
}

}  // namespace smash

#endif  // SRC_INCLUDE_NUMERICS_H_
