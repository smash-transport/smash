/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_IOMANIPULATORS_H_
#define SRC_INCLUDE_IOMANIPULATORS_H_

#include <ostream>

namespace smash {

/**
 * \ingroup logging
 *
 * Stream modifier to align the next object to a specific width \p w.
 *
 * \tparam w The number of characters the field should have in the output.
 * \tparam p The floating precision.
 * \tparam CharT Character type of the output stream.
 * \tparam Traits Traits of the output stream.
 * \param[inout] s The output stream.
 * \return The output stream.
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

}  // namespace smash

#endif  // SRC_INCLUDE_IOMANIPULATORS_H_
