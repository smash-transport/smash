/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_ALGORITHMS_H_
#define SRC_INCLUDE_ALGORITHMS_H_

namespace Smash {

/**
 * Enforces periodic boundaries on the given collection of values.
 *
 * The implementation assumes that the particle is at most one box length
 * away from the boundary to shift it in.
 * This implies for the velocity of the particles:
 * \f[
 * v_i <= \frac{L}{\Delta t}
 * \f]
 *
 * \param begin Iterator pointing to the first value to check.
 * \param end End iterator.
 * \param length The length of the valid interval.
 *
 * \return A tuple of
 *  * the position inside the box (corrected if needed)
 *  * whether a correction was done
 */
template <typename Iterator>
static bool enforce_periodic_boundaries(
    Iterator begin, const Iterator &end,
    typename std::iterator_traits<Iterator>::value_type length) {
  bool had_to_wrap = false;
  for (; begin != end; ++begin) {
    auto &x = *begin;
    if (x < 0) {
      had_to_wrap = true;
      x += length;
    } else if (x >= length) {
      had_to_wrap = true;
      x -= length;
    }
  }
  return had_to_wrap;
}

}  // namespace Smash

#endif  // SRC_INCLUDE_ALGORITHMS_H_
