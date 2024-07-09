/*
 *
 *    Copyright (c) 2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_TRAITS_H_
#define SRC_INCLUDE_SMASH_TRAITS_H_

#include <sstream>
#include <type_traits>

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

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_TRAITS_H_
