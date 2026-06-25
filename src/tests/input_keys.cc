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
  const Key<int> key{
      a_b, 42, {"0.50"}, [](const int& value) { return value > 0; }};
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
  VERIFY(smearing_cutoff.validate(2.0));
  VERIFY(smearing_cutoff.validate(6.5));
  VERIFY(smearing_cutoff.validate(10.0));
  VERIFY(!smearing_cutoff.validate(1.9));
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

TEST(validators_collision_term) {
  const auto& hf_AQM_b_supp = InputKeys::collTerm_HF_AQMbSuppression;
  VERIFY(hf_AQM_b_supp.validate(0.0));
  VERIFY(hf_AQM_b_supp.validate(1.0));
  VERIFY(!hf_AQM_b_supp.validate(-1.0e-12));
  const auto& hf_AQM_c_supp = InputKeys::collTerm_HF_AQMcSuppression;
  VERIFY(hf_AQM_c_supp.validate(0.0));
  VERIFY(hf_AQM_c_supp.validate(1.0));
  VERIFY(!hf_AQM_c_supp.validate(-1.0e-12));
  const auto& add_elastic = InputKeys::collTerm_additionalElasticCrossSection;
  VERIFY(add_elastic.validate(0.0));
  const auto& collision_criterion = InputKeys::collTerm_collisionCriterion;
  VERIFY(collision_criterion.validate(CollisionCriterion::Covariant));
  const auto& xs_scaling = InputKeys::collTerm_crossSectionScaling;
  VERIFY(xs_scaling.validate(1.0));
  VERIFY(!xs_scaling.validate(0.0));
  VERIFY(!xs_scaling.validate(-1.0));
  const auto& elastic_xs = InputKeys::collTerm_elasticCrossSection;
  VERIFY(elastic_xs.validate(-1.0));
  VERIFY(elastic_xs.validate(0.1));
  VERIFY(!elastic_xs.validate(0.0));
  VERIFY(!elastic_xs.validate(-2.0));
  const auto& elastic_NN_cutoff = InputKeys::collTerm_elasticNNCutoffSqrts;
  VERIFY(elastic_NN_cutoff.validate(1.876));
  VERIFY(elastic_NN_cutoff.validate(2.014));
  VERIFY(!elastic_NN_cutoff.validate(1.875));
  VERIFY(!elastic_NN_cutoff.validate(2.015));
  const auto& min_cell_length = InputKeys::collTerm_fixedMinCellLength;
  VERIFY(min_cell_length.validate(1e-6));
  VERIFY(!min_cell_length.validate(0.0));
  const auto& force_decays = InputKeys::collTerm_forceDecaysAtEnd;
  VERIFY(force_decays.validate(true));
  const auto& decay_initial = InputKeys::collTerm_decayInitial;
  VERIFY(decay_initial.validate(true));
  const auto& included_2_to_2 = InputKeys::collTerm_includedTwoToTwo;
  VERIFY(included_2_to_2.validate(ReactionsBitSet{}.set()));
  const auto& include_decays_end = InputKeys::collTerm_includeDecaysAtTheEnd;
  VERIFY(include_decays_end.validate(false));
  const auto& ignore_decay_width = InputKeys::collTerm_ignoreDecayWidthAtTheEnd;
  VERIFY(ignore_decay_width.validate(false));
  const auto& isotropic = InputKeys::collTerm_isotropic;
  VERIFY(isotropic.validate(false));
  const auto& max_xs = InputKeys::collTerm_maximumCrossSection;
  VERIFY(max_xs.validate(666.0));
  VERIFY(!max_xs.validate(0.0));
  const auto& multi_part_reactions = InputKeys::collTerm_multiParticleReactions;
  VERIFY(multi_part_reactions.validate(MultiParticleReactionsBitSet{}.reset()));
  const auto& nnbar_treatment = InputKeys::collTerm_nnbarTreatment;
  VERIFY(nnbar_treatment.validate(NNbarTreatment::Strings));
  const auto& no_collisions = InputKeys::collTerm_noCollisions;
  VERIFY(no_collisions.validate(false));
  const auto& only_warn = InputKeys::collTerm_onlyWarnForHighProbability;
  VERIFY(only_warn.validate(false));
  const auto& pseudoresonance = InputKeys::collTerm_pseudoresonance;
  VERIFY(pseudoresonance.validate(PseudoResonance::LargestFromUnstable));
  const auto& res_lifetime_mod = InputKeys::collTerm_resonanceLifetimeModifier;
  VERIFY(res_lifetime_mod.validate(1.0));
  VERIFY(!res_lifetime_mod.validate(0.0));
  const auto& spinInteractions = InputKeys::collTerm_spinInteractions;
  VERIFY(spinInteractions.validate(SpinInteractionType::Off));
  const auto& strings = InputKeys::collTerm_strings;
  VERIFY(strings.validate(false));
  const auto& strings_w_prob = InputKeys::collTerm_stringsWithProbability;
  VERIFY(strings_w_prob.validate(true));
  const auto& tot_xs_strategy = InputKeys::collTerm_totXsStrategy;
  VERIFY(tot_xs_strategy.validate(TotalCrossSectionStrategy::TopDownMeasured));
  const auto& two_to_one = InputKeys::collTerm_twoToOne;
  VERIFY(two_to_one.validate(true));
  const auto& use_AQM = InputKeys::collTerm_useAQM;
  VERIFY(use_AQM.validate(true));
  const auto& gauss_cutoff = InputKeys::collTerm_pauliBlocking_gaussianCutoff;
  VERIFY(gauss_cutoff.validate(1.0));
  VERIFY(gauss_cutoff.validate(10.0));
  VERIFY(!gauss_cutoff.validate(0.999));
  VERIFY(!gauss_cutoff.validate(10.001));
  const auto& p_av = InputKeys::collTerm_pauliBlocking_momentumAveragingRadius;
  VERIFY(p_av.validate(0.001));
  VERIFY(!p_av.validate(0.0));
  const auto& sp_av = InputKeys::collTerm_pauliBlocking_spatialAveragingRadius;
  VERIFY(sp_av.validate(0.001));
  VERIFY(!sp_av.validate(0.0));
  const auto& KN_offset = InputKeys::collTerm_stringTrans_KNOffset;
  VERIFY(KN_offset.validate(15.15));
  const auto& pipi_offset = InputKeys::collTerm_stringTrans_pipiOffset;
  VERIFY(pipi_offset.validate(1.12));
  const auto& lower = InputKeys::collTerm_stringTrans_lower;
  VERIFY(lower.validate(0.9));
  const auto& range_NN = InputKeys::collTerm_stringTrans_rangeNN;
  VERIFY(range_NN.validate(std::make_pair(1.876, 2.876)));
  VERIFY(!range_NN.validate(std::make_pair(1.87, 4.5)));
  VERIFY(!range_NN.validate(std::make_pair(4.5, 3.5)));
  const auto& range_Npi = InputKeys::collTerm_stringTrans_rangeNpi;
  VERIFY(range_Npi.validate(std::make_pair(1.076, 2.2)));
  VERIFY(!range_Npi.validate(std::make_pair(1.07, 2.2)));
  VERIFY(!range_Npi.validate(std::make_pair(2.2, 1.5)));
  const auto& range_width = InputKeys::collTerm_stringTrans_range_width;
  VERIFY(range_width.validate(1.0));
  const auto& diquark_supp = InputKeys::collTerm_stringParam_diquarkSuppression;
  VERIFY(diquark_supp.validate(0.0));
  VERIFY(diquark_supp.validate(1.0));
  VERIFY(!diquark_supp.validate(-1e-12));
  const auto& form_time_fact = InputKeys::collTerm_stringParam_formTimeFactor;
  VERIFY(form_time_fact.validate(0.1));
  VERIFY(!form_time_fact.validate(0.0));
  const auto& formation_time = InputKeys::collTerm_stringParam_formationTime;
  VERIFY(formation_time.validate(0.1));
  VERIFY(!formation_time.validate(0.0));
  const auto& gluon_beta = InputKeys::collTerm_stringParam_gluonBeta;
  VERIFY(gluon_beta.validate(0.1));
  VERIFY(!gluon_beta.validate(0.0));
  const auto& gluon_p_min = InputKeys::collTerm_stringParam_gluonPMin;
  VERIFY(gluon_p_min.validate(1e-6));
  VERIFY(!gluon_p_min.validate(0.0));
  const auto& md_ft = InputKeys::collTerm_stringParam_mDependentFormationTimes;
  VERIFY(md_ft.validate(false));
  const auto& quark_alpha = InputKeys::collTerm_stringParam_quarkAlpha;
  VERIFY(quark_alpha.validate(0.1));
  VERIFY(!quark_alpha.validate(0.0));
  const auto& quark_beta = InputKeys::collTerm_stringParam_quarkBeta;
  VERIFY(quark_beta.validate(0.1));
  VERIFY(!quark_beta.validate(0.0));
  const auto& popcorn_rate = InputKeys::collTerm_stringParam_popcornRate;
  VERIFY(popcorn_rate.validate(0.0));
  VERIFY(popcorn_rate.validate(1.0));
  VERIFY(!popcorn_rate.validate(-1e-12));
  VERIFY(!popcorn_rate.validate(1.1));
  const auto& power_pf = InputKeys::collTerm_stringParam_powerParticleFormation;
  VERIFY(power_pf.validate(1.0));
  const auto& prob_P_DUU = InputKeys::collTerm_stringParam_probabilityPToDUU;
  VERIFY(prob_P_DUU.validate(1e-6));
  VERIFY(prob_P_DUU.validate(1.0));
  VERIFY(!prob_P_DUU.validate(0.0));
  VERIFY(!prob_P_DUU.validate(1.1));
  const auto& sep_fb = InputKeys::collTerm_stringParam_separateFragmentBaryon;
  VERIFY(sep_fb.validate(true));
  const auto& sigma_perp = InputKeys::collTerm_stringParam_sigmaPerp;
  VERIFY(sigma_perp.validate(1e-6));
  VERIFY(!sigma_perp.validate(0.0));
  const auto& strange_supp = InputKeys::collTerm_stringParam_strangeSuppression;
  VERIFY(strange_supp.validate(0.0));
  VERIFY(strange_supp.validate(1.0));
  VERIFY(!strange_supp.validate(-1e-12));
  VERIFY(!strange_supp.validate(1.1));
  const auto& string_sigma_t = InputKeys::collTerm_stringParam_stringSigmaT;
  VERIFY(string_sigma_t.validate(0.0));
  VERIFY(string_sigma_t.validate(1.0));
  VERIFY(!string_sigma_t.validate(1.0000000001));
  const auto& string_tension = InputKeys::collTerm_stringParam_stringTension;
  VERIFY(string_tension.validate(0.0));
  VERIFY(!string_tension.validate(-1e-12));
  const auto& string_ZA = InputKeys::collTerm_stringParam_stringZA;
  VERIFY(string_ZA.validate(0.0));
  VERIFY(string_ZA.validate(2.0));
  VERIFY(!string_ZA.validate(-1e-12));
  VERIFY(!string_ZA.validate(2.1));
  const auto& string_ZA_lead = InputKeys::collTerm_stringParam_stringZALeading;
  VERIFY(string_ZA_lead.validate(0.0));
  VERIFY(string_ZA_lead.validate(2.0));
  VERIFY(!string_ZA_lead.validate(-1e-12));
  VERIFY(!string_ZA_lead.validate(2.1));
  const auto& string_ZB = InputKeys::collTerm_stringParam_stringZB;
  VERIFY(string_ZB.validate(0.0));
  VERIFY(string_ZB.validate(2.0));
  VERIFY(!string_ZB.validate(-1e-12));
  VERIFY(!string_ZB.validate(2.1));
  const auto& string_ZB_lead = InputKeys::collTerm_stringParam_stringZBLeading;
  VERIFY(string_ZB_lead.validate(0.0));
  VERIFY(string_ZB_lead.validate(2.0));
  VERIFY(!string_ZB_lead.validate(-1e-12));
  VERIFY(!string_ZB_lead.validate(2.1));
  const auto& use_monash_tune = InputKeys::collTerm_stringParam_useMonashTune;
  VERIFY(use_monash_tune.validate(false));
  const auto& hard_string_transition_mode =
      InputKeys::collTerm_hardStringTransition_mode;
  VERIFY(hard_string_transition_mode.validate(
      HardStringTransitionMode::Exponential));
  VERIFY(hard_string_transition_mode.validate(
      HardStringTransitionMode::Custom_Range));
  const auto& hard_string_transition_energy_range =
      InputKeys::collTerm_hardStringTransition_energyRange;
  VERIFY(hard_string_transition_energy_range.validate({10.0, 20.0}));
  VERIFY(hard_string_transition_energy_range.validate({10.0, 200.0}));
  VERIFY(!hard_string_transition_energy_range.validate({9.9, 20.0}));
  VERIFY(!hard_string_transition_energy_range.validate({10.0, 10.0}));
  VERIFY(!hard_string_transition_energy_range.validate({20.0, 10.0}));
  const auto& dileptons_decays = InputKeys::collTerm_dileptons_decays;
  VERIFY(dileptons_decays.validate(false));
  const auto& photons_2_to_2 = InputKeys::collTerm_photons_twoToTwoScatterings;
  VERIFY(photons_2_to_2.validate(false));
  const auto& photons_brems = InputKeys::collTerm_photons_bremsstrahlung;
  VERIFY(photons_brems.validate(false));
  const auto& fractional_phot = InputKeys::collTerm_photons_fractionalPhotons;
  VERIFY(fractional_phot.validate(1));
  VERIFY(!fractional_phot.validate(0));
}

TEST(validators_modi_collider) {
  const auto& e_kin = InputKeys::modi_collider_eKin;
  VERIFY(e_kin.validate(1.0));
  VERIFY(!e_kin.validate(0.0));
  const auto& e_tot = InputKeys::modi_collider_eTot;
  VERIFY(e_tot.validate(0.123));
  VERIFY(!e_tot.validate(0.0));
  const auto& p_lab = InputKeys::modi_collider_pLab;
  VERIFY(p_lab.validate(1.e-12));
  VERIFY(!p_lab.validate(0.0));
  const auto& sqrt_snn = InputKeys::modi_collider_sqrtSNN;
  VERIFY(sqrt_snn.validate(1.e42));
  VERIFY(!sqrt_snn.validate(0.0));
  const auto& calculation_frame = InputKeys::modi_collider_calculationFrame;
  VERIFY(calculation_frame.validate(CalculationFrame::CenterOfMass));
  const auto& coll_in_nucleus = InputKeys::modi_collider_collisionWithinNucleus;
  VERIFY(coll_in_nucleus.validate(true));
  const auto& fermi_motion = InputKeys::modi_collider_fermiMotion;
  VERIFY(fermi_motion.validate(FermiMotion::Frozen));
  const auto& initial_distance = InputKeys::modi_collider_initialDistance;
  VERIFY(initial_distance.validate(1.0));
  VERIFY(!initial_distance.validate(0.0));
  const auto& proj_diff = InputKeys::modi_collider_projectile_diffusiveness;
  VERIFY(proj_diff.validate(0.0));
  VERIFY(proj_diff.validate(1.0));
  VERIFY(!proj_diff.validate(1.0001));
  const auto& targ_diff = InputKeys::modi_collider_target_diffusiveness;
  VERIFY(targ_diff.validate(0.0));
  VERIFY(targ_diff.validate(1.0));
  VERIFY(!targ_diff.validate(-1.e-12));
  const auto& proj_particles = InputKeys::modi_collider_projectile_particles;
  VERIFY(proj_particles.validate(std::map<PdgCode, int>{{0x2212, 42}}));
  VERIFY(!proj_particles.validate(std::map<PdgCode, int>{{0x2212, -2}}));
  const auto& target_particles = InputKeys::modi_collider_target_particles;
  VERIFY(target_particles.validate(std::map<PdgCode, int>{{0x2112, 666}}));
  VERIFY(!target_particles.validate(std::map<PdgCode, int>{{0x2112, 0}}));
  const auto& projectile_radius = InputKeys::modi_collider_projectile_radius;
  VERIFY(projectile_radius.validate(5.0));
  VERIFY(!projectile_radius.validate(0.0));
  const auto& target_radius = InputKeys::modi_collider_target_radius;
  VERIFY(target_radius.validate(5.0));
  VERIFY(!target_radius.validate(0.0));
  const auto& proj_sat = InputKeys::modi_collider_projectile_saturationDensity;
  VERIFY(proj_sat.validate(0.1));
  VERIFY(proj_sat.validate(0.2));
  VERIFY(!proj_sat.validate(0.000009));
  const auto& targ_sat = InputKeys::modi_collider_target_saturationDensity;
  VERIFY(targ_sat.validate(0.1));
  VERIFY(targ_sat.validate(0.2));
  VERIFY(!targ_sat.validate(0.200001));
  const auto& projectile_e_kin = InputKeys::modi_collider_projectile_eKin;
  VERIFY(projectile_e_kin.validate(987.654321));
  VERIFY(!projectile_e_kin.validate(0.0));
  const auto& target_e_kin = InputKeys::modi_collider_target_eKin;
  VERIFY(target_e_kin.validate(123.456789));
  VERIFY(!target_e_kin.validate(0.0));
  const auto& projectile_e_tot = InputKeys::modi_collider_projectile_eTot;
  VERIFY(projectile_e_tot.validate(987.654321));
  VERIFY(!projectile_e_tot.validate(0.0));
  const auto& target_e_tot = InputKeys::modi_collider_target_eTot;
  VERIFY(target_e_tot.validate(123.456789));
  VERIFY(!target_e_tot.validate(0.0));
  const auto& projectile_p_lab = InputKeys::modi_collider_projectile_pLab;
  VERIFY(projectile_p_lab.validate(111.222));
  VERIFY(!projectile_p_lab.validate(0.0));
  const auto& target_p_lab = InputKeys::modi_collider_target_pLab;
  VERIFY(target_p_lab.validate(444.333));
  VERIFY(!target_p_lab.validate(0.0));
  const auto& proj_d = InputKeys::modi_collider_projectile_custom_fileDirectory;
  VERIFY(proj_d.validate("/tmp"));
  VERIFY(!proj_d.validate("/dev/null"));
  VERIFY(!proj_d.validate("/not-existing-folder"));
  const auto& targ_d = InputKeys::modi_collider_target_custom_fileDirectory;
  VERIFY(targ_d.validate("/tmp"));
  VERIFY(!targ_d.validate("/dev/null"));
  VERIFY(!targ_d.validate("/not-existing-folder"));
  const auto& proj_file = InputKeys::modi_collider_projectile_custom_fileName;
  VERIFY(proj_file.validate("nucleus.txt"));
  VERIFY(!proj_file.validate(""));
  VERIFY(!proj_file.validate("."));
  VERIFY(!proj_file.validate(".."));
  VERIFY(!proj_file.validate("foo/nucleus.txt"));
  const auto& targ_file = InputKeys::modi_collider_target_custom_fileName;
  VERIFY(targ_file.validate("nucleus.txt"));
  VERIFY(!targ_file.validate(""));
  VERIFY(!targ_file.validate("."));
  VERIFY(!targ_file.validate(".."));
  VERIFY(!targ_file.validate("foo/nucleus.txt"));
  const auto& pr_d_aut = InputKeys::modi_collider_projectile_deformed_automatic;
  VERIFY(pr_d_aut.validate(false));
  const auto& ta_def_aut = InputKeys::modi_collider_target_deformed_automatic;
  VERIFY(ta_def_aut.validate(true));
  const auto& pr_d_beta2 = InputKeys::modi_collider_projectile_deformed_beta2;
  VERIFY(pr_d_beta2.validate(-1));
  VERIFY(pr_d_beta2.validate(1));
  VERIFY(!pr_d_beta2.validate(-1.000001));
  const auto& ta_d_beta2 = InputKeys::modi_collider_target_deformed_beta2;
  VERIFY(ta_d_beta2.validate(-1));
  VERIFY(ta_d_beta2.validate(1));
  VERIFY(!ta_d_beta2.validate(1.000001));
  const auto& pr_d_beta3 = InputKeys::modi_collider_projectile_deformed_beta3;
  VERIFY(pr_d_beta3.validate(-1));
  VERIFY(pr_d_beta3.validate(1));
  VERIFY(!pr_d_beta3.validate(-1.000001));
  const auto& ta_d_beta3 = InputKeys::modi_collider_target_deformed_beta3;
  VERIFY(ta_d_beta3.validate(-1));
  VERIFY(ta_d_beta3.validate(1));
  VERIFY(!ta_d_beta3.validate(1.000001));
  const auto& pr_d_beta4 = InputKeys::modi_collider_projectile_deformed_beta4;
  VERIFY(pr_d_beta4.validate(-1));
  VERIFY(pr_d_beta4.validate(1));
  VERIFY(!pr_d_beta4.validate(-1.000001));
  const auto& ta_d_beta4 = InputKeys::modi_collider_target_deformed_beta4;
  VERIFY(ta_d_beta4.validate(-1));
  VERIFY(ta_d_beta4.validate(1));
  VERIFY(!ta_d_beta4.validate(1.000001));
  const auto& pr_d_gamma = InputKeys::modi_collider_projectile_deformed_gamma;
  VERIFY(pr_d_gamma.validate(0.0));
  VERIFY(pr_d_gamma.validate(M_PI));
  VERIFY(!pr_d_gamma.validate(3.1416));
  const auto& ta_d_gamma = InputKeys::modi_collider_target_deformed_gamma;
  VERIFY(ta_d_gamma.validate(0.0));
  VERIFY(ta_d_gamma.validate(M_PI));
  VERIFY(!ta_d_gamma.validate(3.1416));
  const auto& projectile_alpha_clustered_automatic =
      InputKeys::modi_collider_projectile_alphaClustered_automatic;
  VERIFY(projectile_alpha_clustered_automatic.validate(true));
  const auto& target_alpha_clustered_automatic =
      InputKeys::modi_collider_target_alphaClustered_automatic;
  VERIFY(target_alpha_clustered_automatic.validate(false));
  const auto& projectile_alpha_side_length =
      InputKeys::modi_collider_projectile_alphaClustered_sideLength;
  VERIFY(projectile_alpha_side_length.validate(1.0));
  VERIFY(!projectile_alpha_side_length.validate(0.0));
  const auto& target_alpha_side_length =
      InputKeys::modi_collider_target_alphaClustered_sideLength;
  VERIFY(target_alpha_side_length.validate(1.0));
  VERIFY(!target_alpha_side_length.validate(0.0));
  const auto& proj_phi = InputKeys::modi_collider_projectile_orientation_phi;
  VERIFY(proj_phi.validate(0.0));
  VERIFY(proj_phi.validate(2 * M_PI));
  VERIFY(!proj_phi.validate(-0.1));
  const auto& targ_phi = InputKeys::modi_collider_target_orientation_phi;
  VERIFY(targ_phi.validate(0.0));
  VERIFY(targ_phi.validate(2 * M_PI));
  VERIFY(!targ_phi.validate(2.000001 * M_PI));
  const auto& pr_theta = InputKeys::modi_collider_projectile_orientation_theta;
  VERIFY(pr_theta.validate(0.0));
  VERIFY(pr_theta.validate(M_PI));
  VERIFY(!pr_theta.validate(-0.1));
  const auto& targ_theta = InputKeys::modi_collider_target_orientation_theta;
  VERIFY(targ_theta.validate(0.0));
  VERIFY(targ_theta.validate(M_PI));
  VERIFY(!targ_theta.validate(-1.e-10 * M_PI));
  const auto& proj_psi = InputKeys::modi_collider_projectile_orientation_psi;
  VERIFY(proj_psi.validate(0.0));
  VERIFY(proj_psi.validate(2 * M_PI));
  VERIFY(!proj_psi.validate(2.000001 * M_PI));
  const auto& targ_psi = InputKeys::modi_collider_target_orientation_psi;
  VERIFY(targ_psi.validate(0.0));
  VERIFY(targ_psi.validate(2 * M_PI));
  VERIFY(!targ_psi.validate(-1.e-11 * M_PI));
  const auto& pr_rnd = InputKeys::modi_collider_projectile_orientation_randRot;
  VERIFY(pr_rnd.validate(true));
  const auto& targ_rnd = InputKeys::modi_collider_target_orientation_randRot;
  VERIFY(targ_rnd.validate(false));
  const auto& impact_max = InputKeys::modi_collider_impact_max;
  VERIFY(impact_max.validate(0.0));
  VERIFY(!impact_max.validate(-1.e-12));
  const auto& imp_rnd_pl = InputKeys::modi_collider_impact_randomReactionPlane;
  VERIFY(imp_rnd_pl.validate(true));
  const auto& impact_range = InputKeys::modi_collider_impact_range;
  VERIFY(impact_range.validate({0.0, 0.0}));
  VERIFY(!impact_range.validate({-1.e-12, 2.0}));
  VERIFY(!impact_range.validate({1.0, -17.0}));
  const auto& impact_sample = InputKeys::modi_collider_impact_sample;
  VERIFY(impact_sample.validate(Sampling::Uniform));
  const auto& impact_value = InputKeys::modi_collider_impact_value;
  VERIFY(impact_value.validate(0.0));
  VERIFY(!impact_value.validate(-1.e-17));
  const auto& impact_values = InputKeys::modi_collider_impact_values;
  VERIFY(impact_values.validate({0.0, 1.0, 2.0}));
  VERIFY(!impact_values.validate({-1.0, 0.0, 1.0}));
  const auto& impact_yields = InputKeys::modi_collider_impact_yields;
  VERIFY(impact_yields.validate({0.0, 1.0, 2.0}));
  VERIFY(!impact_yields.validate({-1.0, 1.0}));
  const auto& ic_type = InputKeys::modi_collider_initialConditions_type;
  VERIFY(ic_type.validate(FluidizationType::ConstantTau));
  const auto& ic_lower = InputKeys::modi_collider_initialConditions_lowerBound;
  VERIFY(ic_lower.validate(1.0));
  VERIFY(!ic_lower.validate(0.0));
  const auto& ic_prop_t = InputKeys::modi_collider_initialConditions_properTime;
  VERIFY(ic_prop_t.validate(1.0));
  VERIFY(!ic_prop_t.validate(0.0));
  const auto& ic_prop_t_sc = InputKeys::modi_collider_initialConditions_scaling;
  VERIFY(ic_prop_t_sc.validate(1.0));
  VERIFY(!ic_prop_t_sc.validate(0.0));
  const auto& ic_pT_cut = InputKeys::modi_collider_initialConditions_pTCut;
  VERIFY(ic_pT_cut.validate(0.0));
  VERIFY(ic_pT_cut.validate(1.0));
  VERIFY(!ic_pT_cut.validate(-0.1));
  const auto& ic_cut = InputKeys::modi_collider_initialConditions_rapidityCut;
  VERIFY(ic_cut.validate(0.0));
  VERIFY(ic_cut.validate(1.0));
  VERIFY(!ic_cut.validate(-0.1));
  const auto& en_thr = InputKeys::modi_collider_initialConditions_eDenThreshold;
  VERIFY(en_thr.validate(0.1));
  VERIFY(!en_thr.validate(0.0));
  const auto& min_time = InputKeys::modi_collider_initialConditions_minTime;
  VERIFY(min_time.validate(0.0));
  VERIFY(!min_time.validate(-1.e-17));
  const auto& max_time = InputKeys::modi_collider_initialConditions_maxTime;
  VERIFY(max_time.validate(1.e6));
  VERIFY(!max_time.validate(0.0));
  const auto& fl_cells = InputKeys::modi_collider_initialConditions_fluidCells;
  VERIFY(fl_cells.validate(2));
  VERIFY(!fl_cells.validate(1));
  const auto& fl_pr = InputKeys::modi_collider_initialConditions_fluidProcesses;
  VERIFY(fl_pr.validate(FluidizableProcessesBitSet{}.set(0)));
  const auto& delay_initial_elastic =
      InputKeys::modi_collider_initialConditions_delayInitialElastic;
  VERIFY(delay_initial_elastic.validate(false));
  const auto& form_time_fraction =
      InputKeys::modi_collider_initialConditions_formTimeFraction;
  VERIFY(form_time_fraction.validate(0.0));
  VERIFY(!form_time_fraction.validate(-1.e-12));
}

TEST(validators_modi_sphere) {
  const auto& init_mult = InputKeys::modi_sphere_initialMultiplicities;
  VERIFY(init_mult.validate(std::map<PdgCode, int>{{{0x2212}, 10}}));
  VERIFY(!init_mult.validate(std::map<PdgCode, int>{{{0x2212}, 0}}));
  VERIFY(!init_mult.validate(std::map<PdgCode, int>{{{0x2212}, -1}}));
  const auto& radius = InputKeys::modi_sphere_radius;
  VERIFY(radius.validate(1.0));
  VERIFY(!radius.validate(0.0));
  VERIFY(!radius.validate(-1.0));
  const auto& start_time = InputKeys::modi_sphere_startTime;
  VERIFY(start_time.validate(0.0));
  VERIFY(start_time.validate(-3.14e66));
  const auto& temperature = InputKeys::modi_sphere_temperature;
  VERIFY(temperature.validate(0.1));
  VERIFY(!temperature.validate(0.0));
  const auto& acc_res_widths = InputKeys::modi_sphere_accountResonanceWidths;
  VERIFY(acc_res_widths.validate(true));
  const auto& add_radial_velocity = InputKeys::modi_sphere_addRadialVelocity;
  VERIFY(add_radial_velocity.validate(0.0));
  VERIFY(add_radial_velocity.validate(1.0));
  VERIFY(!add_radial_velocity.validate(1.1));
  VERIFY(!add_radial_velocity.validate(-1.e-123));
  const auto& add_vel_exp = InputKeys::modi_sphere_addRadialVelocityExponent;
  VERIFY(add_vel_exp.validate(0.0));
  VERIFY(add_vel_exp.validate(2.5));
  VERIFY(!add_vel_exp.validate(-0.1));
  const auto& baryon_chem_pot = InputKeys::modi_sphere_baryonChemicalPotential;
  VERIFY(baryon_chem_pot.validate(-1.0));
  const auto& charge_chem_pot = InputKeys::modi_sphere_chargeChemicalPotential;
  VERIFY(charge_chem_pot.validate(-2.0));
  const auto& ic = InputKeys::modi_sphere_initialCondition;
  VERIFY(ic.validate(SphereInitialCondition::ThermalMomentaBoltzmann));
  const auto& s_chem_pot = InputKeys::modi_sphere_strangeChemicalPotential;
  VERIFY(s_chem_pot.validate(1.e10));
  const auto& hf_multiplier = InputKeys::modi_sphere_heavyFlavorMultiplier;
  VERIFY(hf_multiplier.validate(1.0));
  const auto& use_therm_mult = InputKeys::modi_sphere_useThermalMultiplicities;
  VERIFY(use_therm_mult.validate(true));
  const auto& jet_pdg = InputKeys::modi_sphere_jet_jetPdg;
  VERIFY(jet_pdg.validate(PdgCode{0x2212}));
  const auto& jet_momentum = InputKeys::modi_sphere_jet_jetMomentum;
  VERIFY(jet_momentum.validate(0.1));
  VERIFY(!jet_momentum.validate(0.0));
  const auto& jet_position = InputKeys::modi_sphere_jet_jetPosition;
  VERIFY(jet_position.validate(std::array<double, 3>{{0.0, -1.23, 1.e20}}));
  const auto& jet_back_to_back = InputKeys::modi_sphere_jet_backToBack;
  VERIFY(jet_back_to_back.validate(true));
  const auto& b2b_separation = InputKeys::modi_sphere_jet_backToBackSeparation;
  VERIFY(b2b_separation.validate(0.01));
  VERIFY(!b2b_separation.validate(0.0));
  VERIFY(!b2b_separation.validate(-0.1));
}

TEST(validators_modi_box) {
  const auto& init_mult = InputKeys::modi_box_initialMultiplicities;
  VERIFY(init_mult.validate(std::map<PdgCode, int>{{{0x2212}, 10}}));
  VERIFY(!init_mult.validate(std::map<PdgCode, int>{{{0x2212}, 0}}));
  VERIFY(!init_mult.validate(std::map<PdgCode, int>{{{0x2212}, -1}}));
  const auto& ic = InputKeys::modi_box_initialCondition;
  VERIFY(ic.validate(BoxInitialCondition::ThermalMomentaBoltzmann));
  const auto& length = InputKeys::modi_box_length;
  VERIFY(length.validate(1.0));
  VERIFY(!length.validate(0.0));
  const auto& start_time = InputKeys::modi_box_startTime;
  VERIFY(start_time.validate(0.0));
  VERIFY(start_time.validate(-3.14e66));
  const auto& temperature = InputKeys::modi_box_temperature;
  VERIFY(temperature.validate(0.1));
  VERIFY(!temperature.validate(0.0));
  VERIFY(!temperature.validate(-0.1));
  const auto& acc_res_widths = InputKeys::modi_box_accountResonanceWidths;
  VERIFY(acc_res_widths.validate(true));
  const auto& baryon_chem_pot = InputKeys::modi_box_baryonChemicalPotential;
  VERIFY(baryon_chem_pot.validate(-1.0));
  const auto& charge_chem_pot = InputKeys::modi_box_chargeChemicalPotential;
  VERIFY(charge_chem_pot.validate(-2.0));
  const auto& equilibration_time = InputKeys::modi_box_equilibrationTime;
  VERIFY(equilibration_time.validate(-666.6));
  VERIFY(equilibration_time.validate(0.0));
  VERIFY(equilibration_time.validate(1.5e27));
  const auto& strange_chem_pot = InputKeys::modi_box_strangeChemicalPotential;
  VERIFY(strange_chem_pot.validate(1.e10));
  const auto& use_therm_mult = InputKeys::modi_box_useThermalMultiplicities;
  VERIFY(use_therm_mult.validate(true));
  const auto& jet_momentum = InputKeys::modi_box_jet_jetMomentum;
  VERIFY(jet_momentum.validate(0.1));
  VERIFY(!jet_momentum.validate(0.0));
  const auto& jet_pdg = InputKeys::modi_box_jet_jetPdg;
  VERIFY(jet_pdg.validate(PdgCode{0x2212}));
}

TEST(validators_modi_list) {
  const auto& file_directory = InputKeys::modi_list_fileDirectory;
  VERIFY(file_directory.validate("/tmp"));
  VERIFY(!file_directory.validate("/dev/null"));
  VERIFY(!file_directory.validate("/not-existing-directory"));
  const auto& filename = InputKeys::modi_list_filename;
  VERIFY(filename.validate("dummy.txt"));
  VERIFY(!filename.validate(""));
  VERIFY(!filename.validate("."));
  VERIFY(!filename.validate(".."));
  VERIFY(!filename.validate("foo/bar.dat"));
  const auto& file_prefix = InputKeys::modi_list_filePrefix;
  VERIFY(file_prefix.validate("prefix"));
  const auto& shift_id = InputKeys::modi_list_shiftId;
  VERIFY(shift_id.validate(42));
  VERIFY(shift_id.validate(-17));
  const auto& opt_quantities = InputKeys::modi_list_optionalQuantities;
  VERIFY(opt_quantities.validate(std::vector<std::string>{"ID", "charge"}));
  VERIFY(opt_quantities.validate(std::vector<std::string>{"proc_type"}));
  VERIFY(opt_quantities.validate(std::vector<std::string>{}));
  VERIFY(!opt_quantities.validate(std::vector<std::string>{""}));
}

TEST(validators_modi_listbox) {
  const auto& file_directory = InputKeys::modi_listBox_fileDirectory;
  VERIFY(file_directory.validate("/tmp"));
  VERIFY(!file_directory.validate("/dev/null"));
  VERIFY(!file_directory.validate("/not-existing-directory"));
  const auto& filename = InputKeys::modi_listBox_filename;
  VERIFY(filename.validate("dummy.txt"));
  VERIFY(!filename.validate(""));
  VERIFY(!filename.validate("."));
  VERIFY(!filename.validate(".."));
  VERIFY(!filename.validate("foo/bar.dat"));
  const auto& file_prefix = InputKeys::modi_listBox_filePrefix;
  VERIFY(file_prefix.validate("prefix"));
  const auto& length = InputKeys::modi_listBox_length;
  VERIFY(length.validate(1.0));
  VERIFY(!length.validate(0.0));
  VERIFY(!length.validate(-1.0));
  const auto& shift_id = InputKeys::modi_listBox_shiftId;
  VERIFY(shift_id.validate(0));
  VERIFY(shift_id.validate(100));
  VERIFY(shift_id.validate(-1));
  const auto& opt_quantities = InputKeys::modi_listBox_optionalQuantities;
  VERIFY(opt_quantities.validate(std::vector<std::string>{"ID", "charge"}));
  VERIFY(opt_quantities.validate(std::vector<std::string>{"ncoll"}));
  VERIFY(opt_quantities.validate(std::vector<std::string>{}));
  VERIFY(!opt_quantities.validate(std::vector<std::string>{"ID", "invalid"}));
}

TEST(validators_output) {
  const auto& density_type = InputKeys::output_densityType;
  VERIFY(density_type.validate(DensityType::Baryon));
  const auto& output_interval = InputKeys::output_outputInterval;
  VERIFY(output_interval.validate(1.0));
  VERIFY(!output_interval.validate(0.0));
  const auto& output_times = InputKeys::output_outputTimes;
  VERIFY(output_times.validate({-0.1, 0.0, 1.0, 2.0, 10.0}));
  const auto& particles_format = InputKeys::output_particles_format;
  VERIFY(particles_format.validate({"None"}));
  VERIFY(particles_format.validate({"ASCII"}));
  VERIFY(particles_format.validate({"Binary"}));
  VERIFY(particles_format.validate({"HepMC"}));
  VERIFY(particles_format.validate({"Root"}));
  VERIFY(particles_format.validate({"VTK"}));
  VERIFY(!particles_format.validate({"Lattice_ASCII"}));
  VERIFY(!particles_format.validate({"None", "ASCII"}));
  VERIFY(!particles_format.validate({"InvalidFormat"}));
  const auto& collisions_format = InputKeys::output_collisions_format;
  VERIFY(collisions_format.validate({"None"}));
  VERIFY(collisions_format.validate({"ASCII"}));
  VERIFY(collisions_format.validate({"Binary"}));
  VERIFY(collisions_format.validate({"HepMC"}));
  VERIFY(collisions_format.validate({"Root"}));
  VERIFY(!collisions_format.validate({"VTK"}));
  VERIFY(!collisions_format.validate({"None", "ASCII"}));
  VERIFY(!collisions_format.validate({"For_vHLLE"}));
  const auto& dileptons_format = InputKeys::output_dileptons_format;
  VERIFY(dileptons_format.validate({"None"}));
  VERIFY(dileptons_format.validate({"ASCII"}));
  VERIFY(dileptons_format.validate({"Binary"}));
  VERIFY(dileptons_format.validate({"Root"}));
  VERIFY(!dileptons_format.validate({"HepMC"}));
  VERIFY(!dileptons_format.validate({"VTK"}));
  VERIFY(!dileptons_format.validate({"None", "ASCII"}));
  VERIFY(!dileptons_format.validate({"InvalidFormat"}));
  const auto& photons_format = InputKeys::output_photons_format;
  VERIFY(photons_format.validate({"None"}));
  VERIFY(photons_format.validate({"ASCII"}));
  VERIFY(photons_format.validate({"Binary"}));
  VERIFY(photons_format.validate({"Root"}));
  VERIFY(!photons_format.validate({"HepMC"}));
  VERIFY(!photons_format.validate({"VTK"}));
  VERIFY(!photons_format.validate({"None", "ASCII"}));
  VERIFY(!photons_format.validate({"InvalidFormat"}));
  const auto& ic_format = InputKeys::output_initialConditions_format;
  VERIFY(ic_format.validate({"None"}));
  VERIFY(ic_format.validate({"For_vHLLE"}));
  VERIFY(ic_format.validate({"ASCII"}));
  VERIFY(ic_format.validate({"Binary"}));
  VERIFY(ic_format.validate({"Root"}));
  VERIFY(!ic_format.validate({"HepMC"}));
  VERIFY(!ic_format.validate({"VTK"}));
  VERIFY(!ic_format.validate({"None", "ASCII"}));
  VERIFY(!ic_format.validate({"InvalidFormat"}));
  const auto& rivet_format = InputKeys::output_rivet_format;
  VERIFY(rivet_format.validate({"None"}));
  VERIFY(rivet_format.validate({"YODA"}));
  VERIFY(rivet_format.validate({"YODA-full"}));
  VERIFY(!rivet_format.validate({"ASCII"}));
  VERIFY(!rivet_format.validate({"Binary"}));
  VERIFY(!rivet_format.validate({"HepMC"}));
  VERIFY(!rivet_format.validate({"Root"}));
  VERIFY(!rivet_format.validate({"VTK"}));
  VERIFY(!rivet_format.validate({"None", "YODA"}));
  VERIFY(!rivet_format.validate({"InvalidFormat"}));
  const auto& coulomb_format = InputKeys::output_coulomb_format;
  VERIFY(coulomb_format.validate({"None"}));
  VERIFY(coulomb_format.validate({"VTK"}));
  VERIFY(!coulomb_format.validate({"ASCII"}));
  VERIFY(!coulomb_format.validate({"Binary"}));
  VERIFY(!coulomb_format.validate({"HepMC"}));
  VERIFY(!coulomb_format.validate({"Root"}));
  VERIFY(!coulomb_format.validate({"None", "YODA"}));
  VERIFY(!coulomb_format.validate({"InvalidFormat"}));
  const auto& thermodynamics_format = InputKeys::output_thermodynamics_format;
  VERIFY(thermodynamics_format.validate({"None"}));
  VERIFY(thermodynamics_format.validate({"Lattice_ASCII"}));
  VERIFY(thermodynamics_format.validate({"Lattice_Binary"}));
  VERIFY(thermodynamics_format.validate({"VTK"}));
  VERIFY(!thermodynamics_format.validate({"ASCII"}));
  VERIFY(!thermodynamics_format.validate({"Binary"}));
  VERIFY(!thermodynamics_format.validate({"HepMC"}));
  VERIFY(!thermodynamics_format.validate({"Root"}));
  VERIFY(!thermodynamics_format.validate({"None", "VTK"}));
  VERIFY(!thermodynamics_format.validate({"InvalidFormat"}));
  const auto& particles_extended = InputKeys::output_particles_extended;
  VERIFY(particles_extended.validate(true));
  const auto& particles_quantities = InputKeys::output_particles_quantities;
  VERIFY(particles_quantities.validate({"px", "py", "pz"}));
  VERIFY(!particles_quantities.validate({"px", "invalid_quantity"}));
  const auto& particles_only_final = InputKeys::output_particles_onlyFinal;
  VERIFY(particles_only_final.validate(OutputOnlyFinal::Yes));
  const auto& collisions_extended = InputKeys::output_collisions_extended;
  VERIFY(collisions_extended.validate(true));
  const auto& collisions_quantities = InputKeys::output_collisions_quantities;
  VERIFY(collisions_quantities.validate({"t", "x"}));
  VERIFY(!collisions_quantities.validate({"px", "Foo"}));
  const auto& coll_print_start_end = InputKeys::output_collisions_printStartEnd;
  VERIFY(coll_print_start_end.validate(true));
  const auto& dileptons_extended = InputKeys::output_dileptons_extended;
  VERIFY(dileptons_extended.validate(true));
  const auto& dileptons_quantities = InputKeys::output_dileptons_quantities;
  VERIFY(dileptons_quantities.validate({"mass"}));
  VERIFY(!dileptons_quantities.validate({"px", "Foo"}));
  const auto& photons_extended = InputKeys::output_photons_extended;
  VERIFY(photons_extended.validate(true));
  const auto& photons_quantities = InputKeys::output_photons_quantities;
  VERIFY(photons_quantities.validate({"charge"}));
  VERIFY(!photons_quantities.validate({"px", "Foo"}));
  const auto& ic_extended = InputKeys::output_initialConditions_extended;
  VERIFY(ic_extended.validate(true));
  const auto& ic_quantities = InputKeys::output_initialConditions_quantities;
  VERIFY(ic_quantities.validate({"mt"}));
  VERIFY(!ic_quantities.validate({"px", "Foo"}));
  const auto& ic_lower_bound = InputKeys::output_initialConditions_lowerBound;
  VERIFY(ic_lower_bound.validate(0.5));
  const auto& ic_proper_time = InputKeys::output_initialConditions_properTime;
  VERIFY(ic_proper_time.validate(1.0));
  const auto& ic_pt_cut = InputKeys::output_initialConditions_pTCut;
  VERIFY(ic_pt_cut.validate(0.5));
  const auto& ic_rapidity_cut = InputKeys::output_initialConditions_rapidityCut;
  VERIFY(ic_rapidity_cut.validate(1.0));
  const auto& rivet_analyses = InputKeys::output_rivet_analyses;
  VERIFY(rivet_analyses.validate({"AnalysisName"}));
  const auto& rivet_cross_section = InputKeys::output_rivet_crossSection;
  VERIFY(rivet_cross_section.validate({10.0, 20.0}));
  const auto& rivet_ignore_beams = InputKeys::output_rivet_ignoreBeams;
  VERIFY(rivet_ignore_beams.validate(true));
  const auto& rivet_logging = InputKeys::output_rivet_logging;
  VERIFY(rivet_logging.validate({{"analysis", "INFO"}}));
  const auto& rivet_paths = InputKeys::output_rivet_paths;
  VERIFY(rivet_paths.validate({"/path/to/rivet"}));
  const auto& rivet_preloads = InputKeys::output_rivet_preloads;
  VERIFY(rivet_preloads.validate({"preload.dat"}));
  const auto& rivet_weights_cap = InputKeys::output_rivet_weights_cap;
  VERIFY(rivet_weights_cap.validate(100.0));
  const auto& rivet_weights_deselect = InputKeys::output_rivet_weights_deselect;
  VERIFY(rivet_weights_deselect.validate({"weight1", "weight2"}));
  const auto& rivet_w_nlo_smear = InputKeys::output_rivet_weights_nloSmearing;
  VERIFY(rivet_w_nlo_smear.validate(0.3));
  const auto& rivet_weights_no_multi = InputKeys::output_rivet_weights_noMulti;
  VERIFY(rivet_weights_no_multi.validate(false));
  const auto& rivet_weights_nominal = InputKeys::output_rivet_weights_nominal;
  VERIFY(rivet_weights_nominal.validate("nominal_weight"));
  const auto& rivet_weights_select = InputKeys::output_rivet_weights_select;
  VERIFY(rivet_weights_select.validate({"weightA"}));
  const auto& th_only_part = InputKeys::output_thermodynamics_onlyParticipants;
  VERIFY(th_only_part.validate(false));
  const auto& th_ignore_unf = InputKeys::output_thermodynamics_ignoreUnformed;
  VERIFY(th_ignore_unf.validate(true));
  const auto& th_position = InputKeys::output_thermodynamics_position;
  VERIFY(th_position.validate({1.0, 2.0, 3.0}));
  const auto& th_quantities = InputKeys::output_thermodynamics_quantites;
  VERIFY(th_quantities.validate({ThermodynamicQuantity::Tmn}));
  const auto& th_smearing = InputKeys::output_thermodynamics_smearing;
  VERIFY(th_smearing.validate(true));
  const auto& th_type = InputKeys::output_thermodynamics_type;
  VERIFY(th_type.validate(DensityType::Baryon));
}

TEST(validators_lattice) {
  const auto& automatic = InputKeys::lattice_automatic;
  VERIFY(automatic.validate(true));
  const auto& cell_number = InputKeys::lattice_cellNumber;
  VERIFY(cell_number.validate({200, 250, 300}));
  VERIFY(!cell_number.validate({0, 1, 1}));
  VERIFY(!cell_number.validate({1, 0, 1}));
  VERIFY(!cell_number.validate({1, 1, 0}));
  VERIFY(!cell_number.validate({-1, 1, 1}));
  const auto& origin = InputKeys::lattice_origin;
  VERIFY(origin.validate({0.0, 0.0, 0.0}));
  VERIFY(origin.validate({-1.5, 2.3, 1.e42}));
  const auto& periodic = InputKeys::lattice_periodic;
  VERIFY(periodic.validate(false));
  const auto& pot_affect_thresh = InputKeys::lattice_potentialsAffectThreshold;
  VERIFY(pot_affect_thresh.validate(false));
  const auto& sizes = InputKeys::lattice_sizes;
  VERIFY(sizes.validate({42.0, 199.9, 1.234}));
  VERIFY(!sizes.validate({0.0, 1.0, 1.0}));
  VERIFY(!sizes.validate({1.0, 0.0, 1.0}));
  VERIFY(!sizes.validate({1.0, 1.0, 0.0}));
  VERIFY(!sizes.validate({-1.0, 1.0, 1.0}));
}

TEST(validators_potentials) {
  const auto& use_out = InputKeys::potentials_use_potentials_outside_lattice;
  VERIFY(use_out.validate(false));
  const auto& skyrme_a = InputKeys::potentials_skyrme_skyrmeA;
  VERIFY(skyrme_a.validate(-1.0));
  VERIFY(!skyrme_a.validate(0.0));
  const auto& skyrme_b = InputKeys::potentials_skyrme_skyrmeB;
  VERIFY(skyrme_b.validate(1.0));
  VERIFY(!skyrme_b.validate(0.0));
  const auto& skyrme_tau = InputKeys::potentials_skyrme_skyrmeTau;
  VERIFY(skyrme_tau.validate(0.6667));
  VERIFY(!skyrme_tau.validate(2.0 / 3.0));
  const auto& symmetry_gamma = InputKeys::potentials_symmetry_gamma;
  VERIFY(symmetry_gamma.validate(0.01));
  VERIFY(!symmetry_gamma.validate(0.0));
  const auto& symmetry_s_pot = InputKeys::potentials_symmetry_sPot;
  VERIFY(symmetry_s_pot.validate(30.0));
  VERIFY(symmetry_s_pot.validate(-10.0));
  const auto& vdf_coeffs = InputKeys::potentials_vdf_coeffs;
  VERIFY(vdf_coeffs.validate({1.0, 2.0, 3.0}));
  VERIFY(vdf_coeffs.validate({}));
  const auto& vdf_powers = InputKeys::potentials_vdf_powers;
  VERIFY(vdf_powers.validate({1.0, 2.0, 3.0}));
  VERIFY(vdf_powers.validate({}));
  VERIFY(!vdf_powers.validate({1.0, -2.0, 3.0}));
  const auto& sat_rho_b = InputKeys::potentials_vdf_satRhoB;
  VERIFY(sat_rho_b.validate(0.13));
  VERIFY(sat_rho_b.validate(0.19));
  VERIFY(!sat_rho_b.validate(0.1299));
  VERIFY(!sat_rho_b.validate(0.1901));
  const auto& r_cut = InputKeys::potentials_coulomb_rCut;
  VERIFY(r_cut.validate(0.1));
  VERIFY(!r_cut.validate(0.0));
  const auto& momentum_dep_c = InputKeys::potentials_momentum_dependence_C;
  VERIFY(momentum_dep_c.validate(0.0));
  VERIFY(momentum_dep_c.validate(-1.e123));
  const auto& mom_dep_lambda = InputKeys::potentials_momentum_dependence_Lambda;
  VERIFY(mom_dep_lambda.validate(1.0));
  VERIFY(!mom_dep_lambda.validate(0.0));
}

TEST(validators_forced_thermalization) {
  const auto& cell_number = InputKeys::forcedThermalization_cellNumber;
  VERIFY(cell_number.validate({100, 100, 200}));
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
