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

#include <map>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "stringify.h"

namespace smash {

/// specialize a type for all of the STL containers.
namespace detail {

/**
 * Implementation of the type trait to infer if a type is an STL container. This
 * is the general case for types that are not STL containers.
 *
 * \attention This is at the moment meant to serve the SMASH codebase only and
 * it is on purpose not to add all STL possible containers here. In particular,
 * the goal is to implement the \c is_writable_to_stream trait in a general way
 * since we need it in the Configuration class through the Key::as_yaml method,
 * where the type streamed is totally generic and we need to avoid using the
 * stream \c operator<< on types that do not support it. If SMASH is used as
 * library and another container is needed, a specialization of the
 * <tt>is_stl_container</tt> template can be easily added in smash::detail
 * namespace.
 *
 * \warning The type \c std::array should NOT be added here as STL container.
 * This would lead to an ambiguity in the \c is_writable_to_stream trait and
 * would break compilation. In fact, it is standard to consider \c std::array as
 * tuple-like and so it is done here.
 *
 * \tparam T The type to be tested.
 */
template <typename T>
struct is_stl_container : std::false_type {};

/**
 * Trait specialization for <tt>std::vector</tt>.
 *
 * \tparam Args The STL container template parameters.
 */
template <typename... Args>
struct is_stl_container<std::vector<Args...>> : std::true_type {};

/**
 * Trait specialization for <tt>std::set</tt>.
 *
 * \tparam Args The STL container template parameters.
 */
template <typename... Args>
struct is_stl_container<std::set<Args...>> : std::true_type {};

/**
 * Trait specialization for <tt>std::map</tt>.
 *
 * \tparam Args The STL container template parameters.
 */
template <typename... Args>
struct is_stl_container<std::map<Args...>> : std::true_type {};

}  // namespace detail

//============================================================================//
/**
 * Type trait to infer if a type is an STL container.  Inheriting from
 * <tt>std::bool_constant</tt> will provide a \c value boolean static constant
 * member set depending on whether the type is an STL container.
 *
 * \tparam T The type to be tested.
 */
template <typename T>
struct is_stl_container
    : std::bool_constant<detail::is_stl_container<std::decay_t<T>>::value> {};

/**
 * Helper alias which is common to be defined next to a type trait.
 */
template <typename T>
inline constexpr bool is_stl_container_v = is_stl_container<T>::value;

//============================================================================//
/**
 * Type trait to infer if a type is tuple-like (\c std::pair or \c std::tuple or
 * \c std::array or few others).
 *
 * \tparam T The type to be tested.
 * \tparam Enable Type set to \c void as utility to implement the type trait.
 */
template <typename T, typename Enable = void>
struct is_tuple_like : std::false_type {};

/**
 * Trait specialization for the case when the type is tuple-like. A tuple-like
 * type is recognised trying to access \c std::tuple_size<T>::value which exists
 * only for tuple-like types.
 *
 * \tparam T The type to be tested.
 */
template <typename T>
struct is_tuple_like<T, std::void_t<decltype(std::tuple_size<T>::value)>>
    : std::true_type {};

/**
 * Helper alias which is common to be defined next to a type trait.
 */
template <typename T>
inline constexpr bool is_tuple_like_v = is_tuple_like<T>::value;

//============================================================================//
/**
 * Type trait to infer if a type is map-like (for the moment only \c std::map is
 * considered). \see detail::is_stl_container for the reason why we limit the
 * considered types.
 *
 * \tparam T The type to be tested.
 */
template <typename T>
struct is_map_like : std::false_type {};

/**
 * Trait specialization for the case when the type is <tt>std::map</tt>.
 *
 * \tparam K The map key type.
 * \tparam V The map value type.
 * \tparam Args The map remaining template types.
 */
template <typename K, typename V, typename... Args>
struct is_map_like<std::map<K, V, Args...>> : std::true_type {};

/**
 * Helper alias which is common to be defined next to a type trait.
 */
template <typename T>
inline constexpr bool is_map_like_v = is_map_like<std::decay_t<T>>::value;

//============================================================================//
/**
 * Type trait to infer if a type can be streamed via the \c << operator. This is
 * the general case for not streamable types.
 *
 * \tparam S Type of the stream.
 * \tparam T Type to be tested.
 * \tparam Enable Type set to \c void as utility to implement the type trait.
 */
template <typename S, typename T, typename Enable = void>
struct is_streamable : std::false_type {};

/**
 * Trait specialization for the case when the type is streamable.
 *
 * \tparam S Type of the stream.
 * \tparam T Type to be tested.
 */
template <typename S, typename T>
struct is_streamable<
    S, T, std::void_t<decltype(std::declval<S&>() << std::declval<T>())>>
    : std::true_type {};

/**
 * Helper alias which is common to be defined next to a type trait.
 */
template <typename S, typename T>
inline constexpr bool is_streamable_v = is_streamable<S, T>::value;

//============================================================================//
/**
 * Type trait to infer if a type can be written via the \c << operator, by that
 * not only meaning that an overload exists, but also that possible contained
 * types can be streamed, too. This is the general case for not streamable
 * types.
 *
 * \tparam S Type of the stream.
 * \tparam T Type to be tested.
 * \tparam Enable Type set to \c void as utility to implement the type trait.
 */
template <typename S, typename T, typename Enable = void>
struct is_writable_to_stream : std::false_type {};

/**
 * Trait specialization for the case in which the type is not a container and
 * not tuple-like. Inheriting from \c std::bool_constant will provide a \c value
 * boolean static constant member set depending on whether the type is
 * streamable to the given stream.
 *
 * \note If \c T is not an STL container, it cannot be map-like and there is no
 * need to explicitly exclude map-like types that are treated in a separate
 * specialization.
 *
 * \tparam S Type of the stream.
 * \tparam T Type to be tested.
 */
template <typename S, typename T>
struct is_writable_to_stream<
    S, T, std::enable_if_t<!is_stl_container_v<T> && !is_tuple_like_v<T>>>
    : std::bool_constant<is_streamable_v<S, T>> {};

/**
 * Trait specialization for the case in which the type is an STL container, but
 * not map-like. Inheriting from \c std::bool_constant will provide a \c value
 * boolean static constant member set depending on whether the type is
 * streamable to the given stream and the container value type is also
 * streamable to the given stream.
 *
 * \note As the container value type might be itself a container or something
 * tuple-like, it is necessary here to recurse on \c T::value_type passing it to
 * the type trait itself.
 *
 * \tparam S Type of the stream.
 * \tparam T Type to be tested.
 */
template <typename S, typename T>
struct is_writable_to_stream<
    S, T, std::enable_if_t<is_stl_container_v<T> && !is_map_like_v<T>>>
    : std::bool_constant<
          is_streamable_v<S, T> &&
          is_writable_to_stream<S, typename T::value_type>::value> {};

/**
 * Trait specialization for the case in which the type is map-like. This time
 * we do not inherit from \c std::bool_constant because we need to define a
 * couple of type aliases to identify the map key and value types. Hence we
 * provide ourselves a \c value member.
 *
 * \note As the map key or value types might be themselves containers or
 * tuple-like, it is necessary also here to recurse on the type trait itself.
 *
 * \tparam S Type of the stream.
 * \tparam T Type to be tested.
 */
template <typename S, typename T>
struct is_writable_to_stream<S, T, std::enable_if_t<is_map_like_v<T>>> {
 private:
  /// Type alias for the key type
  using key_type = typename T::key_type;
  /// Type alias for the value type
  using mapped_type = typename T::mapped_type;

 public:
  /// The type trait value which is set testing the private aliases.
  static constexpr bool value = is_streamable_v<S, T> &&
                                is_writable_to_stream<S, key_type>::value &&
                                is_writable_to_stream<S, mapped_type>::value;
};

/**
 * Trait specialization for the case in which the type is tuple-like. This time
 * we do not inherit from \c std::bool_constant because we need an helper
 * function to check all tuple types. Hence we provide ourselves a \c value
 * member.
 *
 * \note As the tuple value types might be themselves containers or something
 * tuple-like, it is necessary also here to recurse on the type trait itself.
 *
 * \tparam S Type of the stream.
 * \tparam T Type to be tested.
 */
template <typename S, typename T>
struct is_writable_to_stream<S, T, std::enable_if_t<is_tuple_like_v<T>>> {
 private:
  /**
   * Helper function to check whether all tuple types are writable to the given
   * stream. This has to be done with a variadic template, as so is
   * <tt>std::tuple</tt>.
   *
   * \tparam I A pack of \c std::size_t types.
   * \return true if all tuple types are writable to the given stream;
   * \return false otherwise.
   */
  template <std::size_t... I>
  static constexpr bool check(std::index_sequence<I...>) {
    return (is_writable_to_stream<S, std::tuple_element_t<I, T>>::value && ...);
  }

 public:
  /// The type trait value which is set using the private method.
  static constexpr bool value =
      is_streamable_v<S, T> &&
      check(std::make_index_sequence<std::tuple_size_v<T>>{});
};

/**
 * Helper alias which is common to be defined next to a type trait.
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
