/*
 *
 *    Copyright (c) 2022,2024-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/input_keys.h"

#include "smash/stringify.h"

using namespace smash;

TEST(section_compile_time) {
  static constexpr InputSections::Section a{"A"};
  static_assert(a.parent == nullptr);
  static_assert(a.name == "A");
  static constexpr InputSections::Section a_b{"B", &a};
  static_assert(a_b.parent == &a);
  static_assert(a_b.name == "B");
}

TEST(section_concatenation) {
  static constexpr InputSections::Section a{"A"};
  static constexpr InputSections::Section a_b = a + "B";
  static_assert(a_b.parent == &a);
  static_assert(a_b.name == "B");
}

TEST(section_conversion) {
  static constexpr InputSections::Section a{"A"};
  static constexpr InputSections::Section a_b = a + "B";
  static constexpr InputSections::Section a_b_c = a_b + "C";
  const KeyLabels expected_labels = {"A", "B", "C"};
  VERIFY(static_cast<KeyLabels>(a_b_c) == expected_labels);
}

TEST(create_key_from_key_path_conversion) {
  static constexpr InputSections::Section a{"A"};
  static constexpr InputSections::Section a_b = a + "B";
  const Key<int> key{a_b, 42, {"0.50"}};
  VERIFY(key.has_same_labels({"A", "B"}));
  VERIFY(key.default_value() == 42);
}

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
 * least one enum/bitset type contained in the InputKeys::key_references_variant
 * which is needed to be able to set Configuration values using enum types. This
 * test is meant to be a support for the developer when adding new enum for a
 * new Key and possibly forget to implement such a conversion.
 *
 * We need few helper traits:
 *  - one to infer if an enum is to be ignored;
 *  - one to identify if a type is an std::bitset;
 *  - one which contains a method to fold over the parameters pack of a variadic
 *    template and test whether the overload exists, but for enum types only.
 */

// Some enums might belong to deleted keys and hence to be ignored
template <typename T>
struct is_enum_to_be_ignored : std::false_type {};

template <>
struct is_enum_to_be_ignored<RestFrameDensityDerivativesMode> : std::true_type {
};

template <typename T>
inline constexpr bool is_enum_to_be_ignored_v = is_enum_to_be_ignored<T>::value;

template <typename T>
struct is_bitset : std::false_type {};

template <std::size_t N>
struct is_bitset<std::bitset<N>> : std::true_type {};

template <typename T>
inline constexpr bool is_bitset_v = is_bitset<T>::value;

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
           if constexpr (is_bitset_v<typename Ts::type::type>) {
             static_assert(has_to_string_v<typename Ts::type::type>,
                           "Missing to_string overload for this bitset type");
           }
         }(),      // <-- immediately invoked lambda
         void()),  // <-- The comma operator discards the lambda's result and
                   // yields nothing, void(), giving each pack element a common
                   // type for the fold; only side-effects inside the lambda are
                   // intended
     ...);         // <-- fold over all Ts
  }
};

TEST(are_all_enum_keys_convertible_to_string) {
  check_to_string_for_enums<InputKeys::key_references_variant>::validate();
}

TEST(validators_particles_and_decaymodes) {
  VERIFY(InputKeys::particles.validate("dummy"));
  VERIFY(InputKeys::decaymodes.validate("dummy"));
}

TEST(validators_general) {
  VERIFY(InputKeys::gen_endTime.validate(0.5));
  VERIFY(!InputKeys::gen_endTime.validate(0));
  VERIFY(!InputKeys::gen_endTime.validate(-6.6));
  VERIFY(InputKeys::gen_nevents.validate(20));
  VERIFY(!InputKeys::gen_nevents.validate(0));
  VERIFY(!InputKeys::gen_nevents.validate(-5));
  VERIFY(InputKeys::gen_randomseed.validate(0));
  VERIFY(InputKeys::gen_randomseed.validate(1234567890123456789LL));
  VERIFY(InputKeys::gen_randomseed.validate(-1));
  VERIFY(InputKeys::gen_minNonEmptyEnsembles_maximumEnsembles.validate(100));
  VERIFY(!InputKeys::gen_minNonEmptyEnsembles_maximumEnsembles.validate(0));
  VERIFY(!InputKeys::gen_minNonEmptyEnsembles_maximumEnsembles.validate(-1));
  VERIFY(InputKeys::gen_minNonEmptyEnsembles_number.validate(1));
  VERIFY(!InputKeys::gen_minNonEmptyEnsembles_number.validate(0));
  VERIFY(!InputKeys::gen_minNonEmptyEnsembles_number.validate(-1));
  VERIFY(InputKeys::gen_deltaTime.validate(1.0));
  VERIFY(!InputKeys::gen_deltaTime.validate(0.0));
  VERIFY(!InputKeys::gen_deltaTime.validate(-1.0));
  VERIFY(InputKeys::gen_smearingDiscreteWeight.validate(0.15));
  VERIFY(InputKeys::gen_smearingDiscreteWeight.validate(0.99));
  VERIFY(!InputKeys::gen_smearingDiscreteWeight.validate(1.0 / 7.0));
  VERIFY(!InputKeys::gen_smearingDiscreteWeight.validate(1.0));
  VERIFY(!InputKeys::gen_smearingDiscreteWeight.validate(1.1));
  VERIFY(InputKeys::gen_ensembles.validate(1));
  VERIFY(InputKeys::gen_ensembles.validate(10));
  VERIFY(!InputKeys::gen_ensembles.validate(0));
  VERIFY(!InputKeys::gen_ensembles.validate(-1));
  VERIFY(InputKeys::gen_expansionRate.validate(-10.));
  VERIFY(InputKeys::gen_expansionRate.validate(0.));
  VERIFY(InputKeys::gen_expansionRate.validate(10.));
  VERIFY(InputKeys::gen_modus.validate("Box"));
  VERIFY(InputKeys::gen_modus.validate("Collider"));
  VERIFY(InputKeys::gen_modus.validate("List"));
  VERIFY(InputKeys::gen_modus.validate("ListBox"));
  VERIFY(InputKeys::gen_modus.validate("Sphere"));
  VERIFY(!InputKeys::gen_modus.validate("Invalid"));
  VERIFY(InputKeys::gen_derivativesMode.validate(DerivativesMode::Off));
  VERIFY(InputKeys::gen_fieldDerivativesMode.validate(
      FieldDerivativesMode::ChainRule));
  VERIFY(InputKeys::gen_smearingGaussCutoffInSigma.validate(3.0));
  VERIFY(InputKeys::gen_smearingGaussCutoffInSigma.validate(6.5));
  VERIFY(InputKeys::gen_smearingGaussCutoffInSigma.validate(10.0));
  VERIFY(!InputKeys::gen_smearingGaussCutoffInSigma.validate(2.9));
  VERIFY(!InputKeys::gen_smearingGaussCutoffInSigma.validate(10.1));
  VERIFY(InputKeys::gen_smearingGaussianSigma.validate(2.99));
  VERIFY(!InputKeys::gen_smearingGaussianSigma.validate(0.1));
  VERIFY(!InputKeys::gen_smearingGaussianSigma.validate(3.0));
  VERIFY(InputKeys::gen_metricType.validate(ExpansionMode::NoExpansion));
  VERIFY(InputKeys::gen_restFrameDensityDerivativeMode.validate(
      RestFrameDensityDerivativesMode::On));
  VERIFY(InputKeys::gen_smearingMode.validate(SmearingMode::Discrete));
  VERIFY(InputKeys::gen_testparticles.validate(42));
  VERIFY(!InputKeys::gen_testparticles.validate(0));
  VERIFY(!InputKeys::gen_testparticles.validate(-1));
  VERIFY(InputKeys::gen_timeStepMode.validate(TimeStepMode::Fixed));
  VERIFY(InputKeys::gen_smearingTriangularRange.validate(2.345));
  VERIFY(!InputKeys::gen_smearingTriangularRange.validate(0.0));
  VERIFY(!InputKeys::gen_smearingTriangularRange.validate(-3.14));
  VERIFY(InputKeys::gen_useGrid.validate(true));
}

#if 0

// The following code is useful to print all keys in the database for debug
// purposes and it is intentionally left as part of the codebase commented out
// for future needs.

#include "smash/input_keys.h"
#include "smash/traits.h"

template <typename T>
std::enable_if_t<is_writable_to_stream_v<decltype(std::cout), T>> print(
    T const& in) {
  std::cout << std::boolalpha << in << '\n';
}

template <typename T>
std::enable_if_t<!is_writable_to_stream_v<decltype(std::cout), T>> print(
    T const& in) {
  std::cout << "NOT PRINTABLE (" << typeid(in).name() << ")\n";
}

TEST(check_list) {
  for (const auto& key : InputKeys::list) {
    try {
      std::visit(
          [](auto&& arg) {
            std::cout << std::setw(70) << std::string(arg.get()) << "  ";
            if (arg.get().has_dependent_default())
              std::cout << "--> dependent default!\n";
            else {
              const auto def = arg.get().default_value();
              print(def);
            }
          },
          key);
    } catch (const std::bad_optional_access&) {
      std::visit([](auto&&) { std::cout << "--> Mandatory\n"; }, key);
    }
  }
}

#endif