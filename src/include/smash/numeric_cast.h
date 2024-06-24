/*
 *
 *    Copyright (c) 2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_NUMERIC_CAST_H_
#define SRC_INCLUDE_SMASH_NUMERIC_CAST_H_

#include <string>
#include <string_view>

namespace smash {

namespace detail {

/**
 * Get type of variable as string in a human-readable way.
 *
 * @tparam T The type to be returned.
 * @return A \c std::string containing the name of the type.
 */
template <typename T>
constexpr auto type_name() {
  std::string_view name, prefix, suffix;
#ifdef __clang__
  name = __PRETTY_FUNCTION__;
  prefix = "auto smash::detail::type_name() [T = ";
  suffix = "]";
#elif defined(__GNUC__)
  name = __PRETTY_FUNCTION__;
  prefix = "constexpr auto smash::detail::type_name() [with T = ";
  suffix = "]";
#elif defined(_MSC_VER)
  name = __FUNCSIG__;
  prefix = "auto __cdecl smash::detail::type_name<";
  suffix = ">(void)";
#else
  name = "UNKNOWN";
  prefix = "";
  suffix = "";
#endif
  name.remove_prefix(prefix.size());
  name.remove_suffix(suffix.size());
  return std::string{name};
}

}  // namespace detail

/**
 * Function template to perform a safe numeric conversion between types.
 *
 * \tparam To Destination type
 * \tparam From Source type
 * \tparam unnamed The last template parameter makes such that this function
 * template participates in function overload only if the destination type \c To
 * is an arithmetic type (typical usage of SFINAE).
 *
 * \param[in] from Input value to be converted
 * \return The input value converted to the new type \c To .
 * \throw std::domain_error If the value to be converted cannot be represented
 * by any value of the destination type.
 *
 * \note This is adapted from the Microsoft implementation of \c narrow function
 * in the C++ <a href="https://github.com/microsoft/GSL/tree/main">Guidelines
 * Support Library</a>, which is released under MIT license.
 */

template <typename To, class From,
          typename std::enable_if_t<std::is_arithmetic_v<To>, bool> = true>
constexpr To numeric_cast(From from) noexcept(false) {
  /* While this is technically undefined behavior in some cases (i.e., if the
  source value is of floating-point type and cannot fit into the destination
  integral type), the resultant behavior is benign on the platforms that we
  target (i.e., no hardware trap representations are hit). */
  const To to = static_cast<To>(from);

  /*
   * NOTE 1: NaN will always throw, since NaN != NaN
   *
   * NOTE 2: The first condition in the if-clause below is not enough because it
   * might happen that the cast back is matching the initial value when casting
   * signed to unsigned numbers or vice-versa because of the "wrapping around"
   * behaviour. See https://stackoverflow.com/a/52863884 for more information.
   */
  constexpr const bool is_different_signedness =
      (std::is_signed_v<To> != std::is_signed_v<From>);
  if (static_cast<From>(to) != from ||
      (is_different_signedness && ((to < To{}) != (from < From{})))) {
    throw std::domain_error("Numeric cast failed converting '" +
                            detail::type_name<From>() + "' to '" +
                            detail::type_name<To>() + "'.");
  }

  return to;
}

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_NUMERIC_CAST_H_
