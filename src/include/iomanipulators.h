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

template <int w = 9, typename CharT, typename Traits>
inline std::basic_ostream<CharT, Traits> &field(
    std::basic_ostream<CharT, Traits> &s) {
  s.put(s.widen(' '));
  s.setf(std::ios_base::fixed, std::ios_base::floatfield);
  s.width(w);
  s.precision(w - 3);
  return s;
}

}  // namespace Smash

#endif  // SRC_INCLUDE_IOMANIPULATORS_H_
