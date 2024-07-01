
/*
 *
 *    Copyright (c) 2023-2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_CXX17COMPAT_H_
#define SRC_INCLUDE_SMASH_CXX17COMPAT_H_

#include <optional>
#include <type_traits>
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

#if __cplusplus < 202002L
/**
 * Definition for remove_cvref type trait, which is in C++20's standard library
 *
 * \see https://en.cppreference.com/w/cpp/types/remove_cvref
 */
template <class T>
struct remove_cvref {
  /// The type with striped properties
  using type = std::remove_cv_t<std::remove_reference_t<T>>;
};
/**
 * Helper alias which is always defined next to a type trait.
 */
template <class T>
using remove_cvref_t = typename remove_cvref<T>::type;
#else
using std::remove_cvref;
using std::remove_cvref_t;
#endif

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_CXX17COMPAT_H_
