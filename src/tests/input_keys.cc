/*
 *
 *    Copyright (c) 2022,2024-2026
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

/* Here below all validation tests are organized by configuration file section.
 * The philosophy is to test a single valid value and as many invalid values as
 * meaningful for each key, in order to have a good coverage of the validators.
 */
TEST(validators_particles_and_decaymodes) {
  VERIFY(InputKeys::particles.validate("dummy"));
  VERIFY(InputKeys::decaymodes.validate("dummy"));
}

TEST(validators_general) {
  const auto& end_time = InputKeys::gen_endTime;
  VERIFY(end_time.validate(0.5));
  VERIFY(!end_time.validate(0));
  VERIFY(!end_time.validate(-6.6));
  const auto& nevents = InputKeys::gen_nevents;
  VERIFY(nevents.validate(20));
  VERIFY(!nevents.validate(0));
  VERIFY(!nevents.validate(-5));
  const auto& randomseed = InputKeys::gen_randomseed;
  VERIFY(randomseed.validate(0));
  VERIFY(randomseed.validate(1234567890123456789LL));
  VERIFY(randomseed.validate(-1));
  const auto& max_ens = InputKeys::gen_minNonEmptyEnsembles_maximumEnsembles;
  VERIFY(max_ens.validate(100));
  VERIFY(!max_ens.validate(0));
  VERIFY(!max_ens.validate(-1));
  const auto& number = InputKeys::gen_minNonEmptyEnsembles_number;
  VERIFY(number.validate(1));
  VERIFY(!number.validate(0));
  VERIFY(!number.validate(-1));
  const auto& delta_time = InputKeys::gen_deltaTime;
  VERIFY(delta_time.validate(1.0));
  VERIFY(!delta_time.validate(0.0));
  VERIFY(!delta_time.validate(-1.0));
  const auto& smearing_discrete_weight = InputKeys::gen_smearingDiscreteWeight;
  VERIFY(smearing_discrete_weight.validate(0.15));
  VERIFY(smearing_discrete_weight.validate(0.99));
  VERIFY(!smearing_discrete_weight.validate(1.0 / 7.0));
  VERIFY(!smearing_discrete_weight.validate(1.0));
  VERIFY(!smearing_discrete_weight.validate(1.1));
  const auto& ensembles = InputKeys::gen_ensembles;
  VERIFY(ensembles.validate(1));
  VERIFY(ensembles.validate(10));
  VERIFY(!ensembles.validate(0));
  VERIFY(!ensembles.validate(-1));
  const auto& expansion_rate = InputKeys::gen_expansionRate;
  VERIFY(expansion_rate.validate(-10.));
  VERIFY(expansion_rate.validate(0.));
  VERIFY(expansion_rate.validate(10.));
  const auto& modus = InputKeys::gen_modus;
  VERIFY(modus.validate("Box"));
  VERIFY(modus.validate("Collider"));
  VERIFY(modus.validate("List"));
  VERIFY(modus.validate("ListBox"));
  VERIFY(modus.validate("Sphere"));
  VERIFY(!modus.validate("Invalid"));
  const auto& derivatives_mode = InputKeys::gen_derivativesMode;
  VERIFY(derivatives_mode.validate(DerivativesMode::Off));
  const auto& field_derivatives_mode = InputKeys::gen_fieldDerivativesMode;
  VERIFY(field_derivatives_mode.validate(FieldDerivativesMode::ChainRule));
  const auto& smearing_cutoff = InputKeys::gen_smearingGaussCutoffInSigma;
  VERIFY(smearing_cutoff.validate(3.0));
  VERIFY(smearing_cutoff.validate(6.5));
  VERIFY(smearing_cutoff.validate(10.0));
  VERIFY(!smearing_cutoff.validate(2.9));
  VERIFY(!smearing_cutoff.validate(10.1));
  const auto& smearing_sigma = InputKeys::gen_smearingGaussianSigma;
  VERIFY(smearing_sigma.validate(2.99));
  VERIFY(!smearing_sigma.validate(0.1));
  VERIFY(!smearing_sigma.validate(3.0));
  const auto& metric_type = InputKeys::gen_metricType;
  VERIFY(metric_type.validate(ExpansionMode::NoExpansion));
  const auto& derivative_mode = InputKeys::gen_restFrameDensityDerivativeMode;
  VERIFY(derivative_mode.validate(RestFrameDensityDerivativesMode::On));
  const auto& smearing_mode = InputKeys::gen_smearingMode;
  VERIFY(smearing_mode.validate(SmearingMode::Discrete));
  const auto& test_particles = InputKeys::gen_testparticles;
  VERIFY(test_particles.validate(42));
  VERIFY(!test_particles.validate(0));
  VERIFY(!test_particles.validate(-1));
  const auto& time_step_mode = InputKeys::gen_timeStepMode;
  VERIFY(time_step_mode.validate(TimeStepMode::Fixed));
  const auto& smearing_range = InputKeys::gen_smearingTriangularRange;
  VERIFY(smearing_range.validate(2.345));
  VERIFY(!smearing_range.validate(0.0));
  VERIFY(!smearing_range.validate(-3.14));
  const auto& use_grid = InputKeys::gen_useGrid;
  VERIFY(use_grid.validate(true));
}

TEST(validate_logging) {
  const auto& log_default = InputKeys::log_default;
  VERIFY(log_default.validate(einhard::ALL));
  VERIFY(log_default.validate(einhard::TRACE));
  VERIFY(log_default.validate(einhard::DEBUG));
  VERIFY(log_default.validate(einhard::INFO));
  VERIFY(log_default.validate(einhard::WARN));
  VERIFY(log_default.validate(einhard::ERROR));
  VERIFY(log_default.validate(einhard::FATAL));
  VERIFY(log_default.validate(einhard::OFF));
  VERIFY(InputKeys::log_box.validate(einhard::ALL));
  VERIFY(InputKeys::log_collider.validate(einhard::DEBUG));
  VERIFY(InputKeys::log_yamlConfiguration.validate(einhard::ALL));
  VERIFY(InputKeys::log_experiment.validate(einhard::INFO));
  VERIFY(InputKeys::log_grandcanThermalizer.validate(einhard::ALL));
  VERIFY(InputKeys::log_initialConditions.validate(einhard::DEBUG));
  VERIFY(InputKeys::log_list.validate(einhard::INFO));
  VERIFY(InputKeys::log_main.validate(einhard::WARN));
  VERIFY(InputKeys::log_output.validate(einhard::DEBUG));
  VERIFY(InputKeys::log_potentials.validate(einhard::TRACE));
  VERIFY(InputKeys::log_rootsolver.validate(einhard::INFO));
  VERIFY(InputKeys::log_sphere.validate(einhard::ALL));
  VERIFY(InputKeys::log_action.validate(einhard::DEBUG));
  VERIFY(InputKeys::log_clock.validate(einhard::INFO));
  VERIFY(InputKeys::log_crossSections.validate(einhard::ALL));
  VERIFY(InputKeys::log_decayModes.validate(einhard::WARN));
  VERIFY(InputKeys::log_density.validate(einhard::DEBUG));
  VERIFY(InputKeys::log_distributions.validate(einhard::TRACE));
  VERIFY(InputKeys::log_findScatter.validate(einhard::ALL));
  VERIFY(InputKeys::log_fpe.validate(einhard::FATAL));
  VERIFY(InputKeys::log_grid.validate(einhard::DEBUG));
  VERIFY(InputKeys::log_hyperSurfaceCrossing.validate(einhard::INFO));
  VERIFY(InputKeys::log_inputParser.validate(einhard::ALL));
  VERIFY(InputKeys::log_lattice.validate(einhard::DEBUG));
  VERIFY(InputKeys::log_nucleus.validate(einhard::WARN));
  VERIFY(InputKeys::log_particleType.validate(einhard::ALL));
  VERIFY(InputKeys::log_pauliBlocking.validate(einhard::INFO));
  VERIFY(InputKeys::log_propagation.validate(einhard::TRACE));
  VERIFY(InputKeys::log_pythia.validate(einhard::DEBUG));
  VERIFY(InputKeys::log_resonances.validate(einhard::ALL));
  VERIFY(InputKeys::log_scatterAction.validate(einhard::INFO));
  VERIFY(InputKeys::log_scatterActionMulti.validate(einhard::DEBUG));
  VERIFY(InputKeys::log_tmn.validate(einhard::TRACE));
}

TEST(validators_version) {
  const auto& version = InputKeys::version;
  VERIFY(version.validate("1.0"));
  VERIFY(version.validate("3.0"));
  VERIFY(version.validate(""));
  VERIFY(version.validate("any string"));
}

TEST(validators_forced_thermalization) {
  const auto& cell_number = InputKeys::forcedThermalization_cellNumber;
  VERIFY(cell_number.validate({300, 300, 300}));
  VERIFY(!cell_number.validate({0, 1, 1}));
  VERIFY(!cell_number.validate({1, 0, 1}));
  VERIFY(!cell_number.validate({1, 1, 0}));
  VERIFY(!cell_number.validate({-1, 1, 1}));
  const auto& critical_edens = InputKeys::forcedThermalization_criticalEDensity;
  VERIFY(critical_edens.validate(2.0));
  VERIFY(!critical_edens.validate(0.0));
  VERIFY(!critical_edens.validate(2.5));
  const auto& start_time = InputKeys::forcedThermalization_startTime;
  VERIFY(start_time.validate(0.1));
  const auto& time_step = InputKeys::forcedThermalization_timestep;
  VERIFY(time_step.validate(4.0));
  VERIFY(!time_step.validate(0.0));
  VERIFY(!time_step.validate(4.1));
  const auto& algorithm = InputKeys::forcedThermalization_algorithm;
  VERIFY(algorithm.validate(ThermalizationAlgorithm::BiasedBF));
  const auto& lattice_sizes = InputKeys::forcedThermalization_latticeSizes;
  VERIFY(lattice_sizes.validate({1.0, 1.0, 1.0}));
  VERIFY(!lattice_sizes.validate({0.0, 1.0, 1.0}));
  VERIFY(!lattice_sizes.validate({1.0, 0.0, 1.0}));
  VERIFY(!lattice_sizes.validate({1.0, 1.0, 0.0}));
  const auto& microcanonical = InputKeys::forcedThermalization_microcanonical;
  VERIFY(microcanonical.validate(false));
}

#if 0

/* The following code is useful to print all keys in the database for debug
 * purposes and it is intentionally left as part of the codebase commented out
 * for future needs. */

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
