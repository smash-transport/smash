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

#include <algorithm>
#include <cmath>

/**
 * \file
 *
 * Generic algorithms on containers and ranges.
 *
 * This file collects generic algorithms that follow the general idea of C++
 * algorithms as defined in the C++ standard library. These typically work with
 * iterators from arbitrary containers.
 *
 * The C++ standard itself categorizes algorithms into the following:
 * * Non-modifying sequence operations
 * * Mutating sequence operations
 * * Sorting and related operations
 */

namespace Smash {

/**
 * Enforces periodic boundaries on the given collection of values.
 *
 * The values in an arbitrary container, starting from \p begin and ending at \p
 * end, will be checked. If the value is less than 0, \p length will be added to
 * it. If the value is greater than or equal to \p length, \p length will be
 * subtracted from it.
 *
 * The implementation therefore assumes that the values are at most one \p
 * length away from the 0 to \p length range.
 *
 * \param begin Iterator pointing to the first value to check.
 * \param end End iterator.
 * \param length The length of the valid interval.
 *
 * \return whether a correction was done
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

/**
 * Convenience wrapper for \c std::all_of that operates on a complete container.
 *
 * Returns whether all elements in \p c return \c true when passed to \p p.
 * \param c A container of elements to examine.
 * \param p unary predicate.
 */
template <typename Container, typename UnaryPredicate>
inline bool all_of(Container &&c, UnaryPredicate &&p) {
  return std::all_of(std::begin(c), std::end(c),
                     std::forward<UnaryPredicate>(p));
}

}  // namespace Smash

#endif  // SRC_INCLUDE_ALGORITHMS_H_
