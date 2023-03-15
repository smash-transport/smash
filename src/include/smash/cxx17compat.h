
/*
 *
 *    Copyright (c) 2023
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_CXX17COMPAT_H_
#define SRC_INCLUDE_SMASH_CXX17COMPAT_H_

#include <optional>
#include <utility>

namespace smash {

// GCC_COMPILER macro is defined in CMake code when appropriate
#if defined(GCC_COMPILER) && (__GNUC__ == 8 && __GNUC_MINOR__ < 4)
/**
 * An utility function to circumvent GNU compiler bug in some versions.
 *
 * Explicitly specifying the template parameter when using
 * \c std::make_optional is not correctly implemented in GNU compiler
 * from version 8.1.x to version 8.3.x (both included). See
 * https://stackoverflow.com/q/75653199/14967071 for an example.
 *
 * @tparam T The optional type
 * @param value The value to be used to construct the \c std::optional
 * @return A \c std::optional<T> holding \c value
 */
template <typename T>
std::optional<T> make_optional(T&& value) {
  return std::optional<std::decay_t<T>>(std::forward<T>(value));
}
#else
using std::make_optional;
#endif
}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_CXX17COMPAT_H_
