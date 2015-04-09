/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_IOMANIPULATORS_H_
#define SRC_INCLUDE_IOMANIPULATORS_H_

#include <ostream>

namespace Smash {

/**
 * \ingroup logging
 *
 * Stream modifier to align the next object to a specific width \p w.
 *
 * \tparam w The number of characters the field should have in the output.
 */
template <int w = 9, int p = w - 3, typename CharT, typename Traits>
inline std::basic_ostream<CharT, Traits> &field(
    std::basic_ostream<CharT, Traits> &s) {
  s.put(s.widen(' '));
  s.setf(std::ios_base::fixed, std::ios_base::floatfield);
  s.width(w);
  s.precision(p);
  return s;
}

}  // namespace Smash

#endif  // SRC_INCLUDE_IOMANIPULATORS_H_
