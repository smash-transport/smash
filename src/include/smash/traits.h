/*
 *
 *    Copyright (c) 2024-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_TRAITS_H_
#define SRC_INCLUDE_SMASH_TRAITS_H_

#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "stringify.h"

namespace smash {

/**
 * Type trait to infer if a type can be streamed via the \c << operator. This is
 * the general case for not streamable types.
 *
 * @tparam S Type of the stream
 * @tparam T Type to be tested
 * @tparam typename set to void as utility to implement the type trait
 */
template <typename S, typename T, typename = void>
struct is_writable_to_stream : std::false_type {};

/**
 * Trait specialization for the case when the type is streamable.
 *
 * @tparam S Type of the stream
 * @tparam T Type to be tested
 */
template <typename S, typename T>
struct is_writable_to_stream<
    S, T, std::void_t<decltype(std::declval<S&>() << std::declval<T>())>>
    : std::true_type {};

/**
 * Helper alias which is always defined next to a type trait.
 */
template <typename S, typename T>
inline constexpr bool is_writable_to_stream_v =
    is_writable_to_stream<S, T>::value;

/**
 * Type trait to infer if there is an <tt>std::string to_string(T)</tt> overload
 * for a given type <tt>T</tt>. This is the general case in which the overload
 * is missing.
 *
 * \tparam T The type for which the overload has to be checked.
 * \tparam Enable Type set to \c void as utility to implement the type trait.
 */
template <typename T, typename Enable = void>
struct has_to_string : std::false_type {};

/**
 * Trait specialization for the case when the overload is present. The test is
 * done in two parts. The second template parameter is used to test if there
 * exists a \c smash::to_string(T) overload, while inheriting from \c is_same
 * will provide a \c value boolean static constant member set to \c true if the
 * return value is an \c std::string and to \c false otherwise.
 *
 * \attention Here we test the existence of the overload in the \c smash
 * namespace and hence argument-dependent lookup (ADL) does not kick in. Said
 * differently \c smash::to_string is unqualified and non-dependent and
 * non-dependent names are looked up immediately, at the point of the template
 * definition. Therefore, we need to include here the file that has the defined
 * conversions, otherwise compilation would fail.
 *
 * \tparam T The type for which the overload has to be checked.
 */
template <typename T>
struct has_to_string<T,
                     std::void_t<decltype(smash::to_string(std::declval<T>()))>>
    : std::is_same<decltype(smash::to_string(std::declval<T>())), std::string> {
};

/**
 * Trait specialization for \c std::bitset types for which a different signature
 * of the overload is required. Because of how these are used in SMASH when
 * parsing the input YAML file, it makes sense that the \c to_string overload
 * returns an \c std::vector of strings which are the corresponding enum entries
 * converted to string.
 *
 * \tparam N The size of the bitset.
 */
template <std::size_t N>
struct has_to_string<std::bitset<N>, std::void_t<decltype(smash::to_string(
                                         std::declval<std::bitset<N>>()))>>
    : std::is_same<decltype(smash::to_string(std::declval<std::bitset<N>>())),
                   std::vector<std::string>> {};

/**
 * Helper alias which is always defined next to a type trait.
 */
template <typename T>
inline constexpr bool has_to_string_v = has_to_string<T>::value;

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_TRAITS_H_
