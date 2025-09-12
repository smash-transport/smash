/*
 *
 *    Copyright (c) 2022,2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/input_keys.h"

#include "smash/stringify.h"

using namespace smash;

TEST(get_logging_key) {
  auto key = InputKeys::get_logging_key("Main");
  VERIFY(key.has_same_labels({"Logging", "Main"}));
}

TEST_CATCH(get_wrong_logging_key, std::invalid_argument) {
  auto key = InputKeys::get_logging_key("XXX");
}

TEST(get_output_format_key) {
  auto key = InputKeys::get_output_format_key("Particles");
  VERIFY(key.has_same_labels({"Output", "Particles", "Format"}));
}

TEST_CATCH(get_wrong_output_format_key, std::invalid_argument) {
  auto key = InputKeys::get_output_format_key("YYY");
}

/*
 * In the following we want to implement a TEST case that does not do anything,
 * but its compilation should fail if there is no to_string() overload for at
 * least one enum type contained in the InputKeys::key_references_variant which
 * is needed to be able to set Configuration values using enum types. This test
 * is meant to be a support for the developer when adding new enum for a new Key
 * and possibly forget to implement such a conversion.
 *
 * We need an helper trait which contains a method to fold over the parameters
 * pack of a variadic template and test whether the overload exists, but for
 * enum types only.
 */

// Some enums might belong to deleted keys and hence to be ignored
template <typename T>
struct is_enum_to_be_ignored : std::false_type {};

template <>
struct is_enum_to_be_ignored<RestFrameDensityDerivativesMode> : std::true_type {
};

template <typename T>
inline constexpr bool is_enum_to_be_ignored_v = is_enum_to_be_ignored<T>::value;

/* Implement check on all enum types contained in the variant
 *
 * NOTE: Here we do want to be specific and exclusively serve our type, i.e.
 * InputKeys::key_references_variant and not a generic variant. Hence we assume
 * that the types in the variant are std::reference_wrapper containing Key
 * types. Hence, since Ts is the type contained in the std::variant, Ts::type is
 * the type contained in the std::reference_wrapper and Ts::type::type is the
 * type contained in the Key.
 */
template <typename Variant>
struct check_to_string_for_enums;

template <typename... Ts>
struct check_to_string_for_enums<std::variant<Ts...>> {
  static void validate() {
    ((
         [] {
           if constexpr (std::is_enum_v<typename Ts::type::type> &&
                         !is_enum_to_be_ignored_v<typename Ts::type::type>) {
             static_assert(has_to_string_v<typename Ts::type::type>,
                           "Missing to_string overload for this enum type");
           }
         }(),  // <-- immediately invoked lambda
         0),   // <-- turns each lambda call into an int for the fold.
     ...);     // <-- fold over all Ts
  }
};

TEST(are_all_enum_keys_convertible_to_string) {
  check_to_string_for_enums<InputKeys::key_references_variant>::validate();
}
