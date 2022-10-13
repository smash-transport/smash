/*
 *
 *    Copyright (c) 2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_VALIDATION_H_
#define SRC_INCLUDE_SMASH_VALIDATION_H_

#include <any>
#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include "smash/configuration.h"
#include "smash/stringfunctions.h"

namespace smash {
/**
 * Descriptive alias for storing SMASH versions associated to keys metadata.
 * At the moment simply a \c std::string .
 */
using Version = std::string;

/**
 * @brief Object to store a YAML input file key together with metadata
 * associated to it.
 *
 * @note The class is designed such that all keys can be marked as deprecated
 *       and as removed. However, it is not possible to mark a key as removed
 *       without having deprecated it before. A workaround is to deprecate and
 *       remove it in the same version, i.e. specifying the same version twice
 *       at construction.
 *
 * @tparam default_type Type of the key value.
 */
template <typename default_type>
class Key {
 public:
  /**
   * \ingroup exception
   * Thrown when too few or too many versions are passed to the constructor.
   */
  struct WrongNumberOfVersions : public std::runtime_error {
    using std::runtime_error::runtime_error;
  };

  /**
   * @brief Construct a new \c Key object without default value.
   *
   * @param[in] labels The label(s) identifying the key in the YAML input file.
   * @param[in] versions A list of one, two or three version numbers identifying
   * the versions in which the key has been introduced, deprecated and removed,
   * respectively.
   */
  explicit Key(const std::initializer_list<std::string_view>& labels,
               const std::initializer_list<std::string_view>& versions)
      : Key{labels, std::nullopt, versions} {}

  /**
   * @brief Construct a new \c Key object without default value.
   *
   * Note that the default value could be simply taken as parameter of type
   * \c default_type . However, this would complicate delegating construction
   * in the other constructor and the caller can still pass a variable of type
   * \c default_type and the \c std::optional will be constructed without
   * problems.
   *
   * @param[in] labels The label(s) identifying the key in the YAML input file.
   * @param[in] value The key default value.
   * @param[in] versions A list of one, two or three version numbers identifying
   * the versions in which the key has been introduced, deprecated and removed,
   * respectively.
   *
   * @throw WrongNumberOfVersions If \c versions has the wrong size.
   */
  Key(const std::initializer_list<std::string_view>& labels,
      const std::optional<default_type>& value,
      const std::initializer_list<std::string_view>& versions)
      : default_{value}, labels_{labels.begin(), labels.end()} {
    /*
     * The following switch statement is a compact way to initialize the
     * three version member variables without repetition and lot's of logic
     * clauses. The versions variable can have 1, 2 or 3 entries. The use of
     * the iterator is needed, since std::initializer_list has no access
     * operator.
     */
    switch (auto it = versions.end(); versions.size()) {
      case 3:
        removed_in_ = *(--it);
        [[fallthrough]];
      case 2:
        deprecated_in_ = *(--it);
        [[fallthrough]];
      case 1:
        introduced_in_ = *(--it);
        break;
      default:
        throw WrongNumberOfVersions(
            "Key constructor needs one, two or three version numbers.");
    }
  }

  /**
   * @brief Get the default value of the key.
   *
   * @return A \c default_type variable.
   *
   * @throw std::bad_optional_access If the key has no default value.
   */
  default_type default_value() const { return default_.value(); }

  /**
   * @brief Get the SMASH version in which the key has been introduced.
   *
   * @return A \c Version variable.
   */
  Version introduced_in() const noexcept { return introduced_in_; }

  /**
   * @brief Get the SMASH version in which the key has been deprecated.
   *
   * @return A \c Version variable.
   *
   * @throw std::bad_optional_access If the key is not deprecated.
   */
  Version deprecated_in() const { return deprecated_in_.value(); }

  /**
   * @brief Get the SMASH version in which the key has been removed.
   *
   * @return A \c Version variable.
   *
   * @throw std::bad_optional_access If the key is still allowed.
   */
  Version removed_in() const { return removed_in_.value(); }

  /**
   * @brief Get whether the key is deprecated or not.
   *
   * @return \c true if the key is deprecated, \c false otherwise.
   */
  bool is_deprecated() const noexcept { return deprecated_in_.has_value(); }

  /**
   * @brief Get whether the key is still allowed or not.
   *
   * @return \c true if the key is allowed, \c false otherwise.
   */
  bool is_allowed() const noexcept { return !removed_in_.has_value(); }

  /**
   * @brief Check if given labels are the same as those of this object.
   *
   * @param[in] labels Given labels to be checked against.
   *
   * @return \c true if all labels match in the given order,
   * @return \c false otherwise.
   */
  bool has_same_labels(const std::vector<std::string>& labels) const noexcept {
    return std::equal(std::begin(labels_), std::end(labels_),
                      std::begin(labels), std::end(labels));
  }

  /**
   * @brief Converts a Key to a \c std::string using all labels.
   *
   * @return \c std::string with labels concatenated with \c :␣ (colon-space)
   *         and quotes all around.
   */
  explicit operator std::string() const noexcept {
    return smash::quote(smash::join(labels_, ": "));
  }

 private:
  /// SMASH version in which the key has been introduced
  Version introduced_in_{};
  /// SMASH version in which the key has been deprecated, if any
  std::optional<Version> deprecated_in_{};
  /// SMASH version in which the key has been removed, if any
  std::optional<Version> removed_in_{};
  /// Key default value, if any
  std::optional<default_type> default_{};
  /// The label(s) identifying the key in the YAML input file
  std::vector<std::string> labels_{};
};

/*!\Userguide
 * \page input Input
 *
 * There are three input files used by SMASH:
 *
 * - `config.yaml` for configuring the simulation. This file is required. See
 *   \subpage inputconfig.
 * - `particles.txt` for defining the particles used by SMASH. This file is
 *   optional. See \subpage inputparticles.
 * - `decaymodes.txt` for defining the decays (and corresponding resonance
 *   formations) possible in SMASH. This file is
 *   optional. See \subpage inputdecaymodes.
 *
 * \page inputconfig Configuration
 *
 * SMASH is configured via an input file in YAML format. Typically you will
 * start from the supplied `config.yaml` file and modify it according to your
 * needs. If you ever make a mistake there and specify a configuration key that
 * SMASH does not recognize, then on startup it will tell you about the keys it
 * could not make any sense of.
 *
 * By default, SMASH copies the `config.yaml` file used to set up the SMASH run
 * to the output directory of the simulation. For the sake of reproducibility,
 * the randomly generated number seed (if the user specified a negative seed) is
 * inserted into the copied file and the used particles and decay modes are
 * appended as well.
 *
 * \par The available keys are documented on the following pages:
 * \li \subpage input_general_
 * \li \subpage input_logging_
 * \li \subpage input_collision_term_
 * \li \subpage input_modi_
 * \li \subpage input_output_options_
 * \li \subpage input_lattice_
 * \li \subpage input_potentials_
 * \li \subpage input_forced_thermalization_
 *
 * \par Information on formatting of the input file can be found here:
 * \li \subpage input_indentation_
 *
 * \ifnot user
 * \par The relevant functions and classes for input are:
 * \li \ref Configuration
 * \li \ref ExperimentBase::create()
 * \li \ref ColliderModus
 * \li \ref BoxModus
 * \li \ref SphereModus
 * \li \ref ListModus
 * \li \ref ListBoxModus
 * \endif
 */

/*!\Userguide
 * \page configuration_keys Input short reference
 *
 * This is a look-up reference of input keys. Refer to each corresponding page
 * for a detailed description of each key.
 */

/*!\Userguide
 * \page input_general_ General
 *
 * This section in the `config.yaml` file contains all general/global
 * configuration options to SMASH. Before describing all possible keys in
 * detail, let's tart off with a couple of examples.
 *
 * The `General` section in SMASH input file might read as follows:
 *
 *\verbatim
   General:
       Modus: "Collider"
       Delta_Time: 0.1
       Testparticles: 1
       Gaussian_Sigma: 1.0
       Gauss_Cutoff_In_Sigma: 3.0
       End_Time: 100.0
       Randomseed: -1
       Nevents: 20
       Use_Grid: true
       Time_Step_Mode: "Fixed"
   \endverbatim
 *
 * In the case of an expanding sphere setup, change the \key Modus and provide
 * further information about the expansion.
 *\verbatim
      Modus: "Sphere"
      MetricType: "MasslessFRW"
      Expansion_Rate: 0.1
  \endverbatim
 */

/*!\Userguide
 * \page minimum_nonempty_ensembles_ Minimum non-empty ensembles
 *
 * Instead of defining the number of Events it is possible to define a minimum
 * number of ensembles in which an interaction took place. Using this option
 * by providing a `Minimum_Nonempty_Ensembles` section in the input file,
 * events will be calculated until the desired number of non-empty ensembles
 * is generated. If the <tt>\ref nevents_ "Nevents"</tt> key is not specified,
 * <b>this section with all its required keys must be present in the SMASH
 * input file</b>.
 *
 * number of ensembles is equal to the number of events, so that this option
 * will provide the desired number of non-empty events.
 */

/*!\Userguide
 * \page input_logging_ Logging
 *
 * The `Logging` section in the input file controls the logging levels for
 * different areas of the code, each of which can have a different verbosity
 * level. All keys and hence the section itself are optional. Valid key values
 * are the following:
 * - `"ALL"`   &rarr; Log all messages (default)
 * - `"TRACE"` &rarr; The lowest level for messages describing the program flow
 * - `"DEBUG"` &rarr; Debug messages
 * - `"INFO"`  &rarr; Messages of informational nature
 * - `"WARN"`  &rarr; Warning messages
 * - `"ERROR"` &rarr; Non-fatal errors
 * - `"FATAL"` &rarr; Messages that indicate terminal application failure
 * - `"OFF"`   &rarr; If selected no messages will be printed to the output
 *
 * Note that the logging levels `TRACE` and `DEBUG` are only available in
 * debug builds (i.e. running `cmake` with `-DCMAKE_BUILD_TYPE=Debug`).
 */

/**
 * @brief A container to keep track of all ever existed input keys.
 *
 * @remark Each input key exists as static constant member and a reference to it
 *         is stored in the InputKeys::list container. Therefore, the following
 *         steps are needed in order to add a new key.
 *         -# Add a new member being consistent with the existing notation.
 *            Use \c _ to separate YAML sections in the variable name and use
 *            a name that reflects sections. A double underscore in C++ is
 *            reserved and should not be used in identifiers; hence it must not
 *            be used to separate sections. If any label consists of more than
 *            one word, use lowerCamelCase convention, although this violates
 *            the general codebase rules (it adds readability in this case).
 *            Abbreviations are allowed, but be consistent if any already
 *            exists.
 *         -# If the newly introduced key has a new type w.r.t. all existing
 *            keys, you need to add it to the \c key_references_variant alias.
 *            In particular, you need to add a type to the \c std::variant which
 *            will be <tt>std::reference_wrapper<const Key<NEW_TYPE>></tt> with
 *            \c NEW_TYPE replaced by the type of your new key.
 *         -# Add a reference to the newly introduced variable to the
 *            InputKeys::list container. This must be done using \c std::cref as
 *            for the other references. Respecting the members order is welcome.
 *
 * @attention If you need to deprecate or to mark a key as not valid anymore,
 *            add the corresponding SMASH version to the \c Key member
 *            constructor invocation. <b>Do not remove a member if that key
 *            is not valid any more!</b> It is intended to track here keys
 *            that were existing in the past and are not accepted anymore.
 *            Instead, after having added the version in which the key
 *            has been deprecated or removed, <b>adjust the user documentation
 *            by marking the key as deprecated or by removing the key and
 *            its description</b>.
 *
 * @note Ordering of members in this class is imposed by how keys shall appear
 *       in the documentation. In particular, all mandatory keys per page are
 *       listed first and all optional after in respective sub-sections. Every
 *       block of keys contains them <b>in alphabetical order</b>, keep it so.
 *       Although not strictly necessary, all keys belonging to the same page
 *       are put next to each other.
 */
struct InputKeys {
  /*!\Userguide
   * \page input_general_
   * <hr>
   * \par MANDATORY KEYS
   */

  /*!\Userguide
   * \page input_general_
   * \required_key_no_line{end_time_,End_Time,double}
   *
   * The time in fm after which the evolution is stopped. Note
   * that the starting time depends on the chosen `Modus`.
   */
  /**
   * \see_key{end_time_}
   */
  inline static const Key<double> gen_endTime{{"General", "End_Time"}, {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \required_key{modus_,Modus,string}
   *
   * Selects a modus for the calculation, e.g.\ infinite matter
   * calculation, collision of two particles or collision of nuclei. The modus
   * will be configured in \ref input_modi_. Recognized values are:
   * - `"Collider"` &rarr; For collisions of nuclei or compound objects. See
   *   \ref \ColliderModus
   * - `"Sphere"` &rarr; For calculations of the expansion of a thermalized
   *   sphere. See \ref \SphereModus
   * - `"Box"` &rarr; For infinite matter calculation in a rectangular box. See
   *   \ref \BoxModus
   * - `"List"` &rarr; For given external particle list. See \ref \ListModus
   * - `"ListBox"` &rarr; For given external particle list in the Box.
   */
  /**
   * \see_key{modus_}
   */
  inline static const Key<std::string> gen_modus{{"General", "Modus"}, {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \required_key{nevents_,Nevents,int}
   *
   * Number of events to calculate.
   *
   * This key may be omitted on constraint that a minimum number
   * of ensembles containing interactions is requested, see
   * \subpage minimum_nonempty_ensembles_.
   */
  /**
   * \see_key{nevents_}
   */
  inline static const Key<int> gen_nevents{{"General", "Nevents"}, {"1.0"}};

  /*!\Userguide
   * \page minimum_nonempty_ensembles_
   * \required_key{mnee_number_,Number,int}
   *
   * The number of desired non-empty ensembles.\n
   */
  /**
   * \see_key{mnee_number_}
   */
  inline static const Key<int> gen_minNonEmptyEnsembles_number{
      {"General", "Minimum_Nonempty_Ensembles", "Number"}, {"1.0"}};

  /*!\Userguide
   * \page minimum_nonempty_ensembles_
   * \required_key{mnee_maximum_ensembles_,Maximum_Ensembles_Run,int}
   *
   * Maximum number of ensembles run. This number serves as a safeguard
   * against SMASH unexpectedly running for a long time.
   */
  /**
   * \see_key{mnee_maximum_ensembles_}
   */
  inline static const Key<int> gen_minNonEmptyEnsembles_maximumEnsembles{
      {"General", "Minimum_Nonempty_Ensembles", "Maximum_Ensembles_Run"},
      {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \required_key{randomseed_,Randomseed,int}
   *
   * Initial seed for the random number generator. If this is negative, the
   * seed will be randomly generated by the operating system.
   */
  /**
   * \see_key{randomseed_}
   */
  inline static const Key<int> gen_randomseed{{"General", "Randomseed"},
                                              {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * <hr>
   * \par OPTIONAL KEYS
   */

  /*!\Userguide
   * \page input_general_
   * \optional_key_no_line{delta_time_,Delta_Time,double,1.0}
   *
   * Fixed time step at which the collision-finding grid is recreated, and, if
   * potentials are on, momenta are updated according to the equations of
   * motion. The collision-finding grid finds all the collisions from time
   * t_{beginning_of_timestep} until time t_{beginning_of_timestep} +
   * Delta_Time, and puts them into a vector. The collisions are then sorted in
   * order of occurrence, and particles are propagated from collision to
   * collision. After each performed collision, additional collisions are found
   * for outgoing particles and merged into the sorted vector.
   *
   * If potentials are on, the Delta_Time should be small enough, typically
   * around 0.1 fm/c. However, if potentials are off, it can be arbitrarily
   * large. In this case it only influences the runtime, but not physics.
   * If Time_Step_Mode = None is chosen, then the user-provided value of
   * Delta_Time is ignored and Delta_Time is set to the End_Time.
   */
  /**
   * \see_key{delta_time_}
   */
  inline static const Key<double> gen_deltaTime{
      {"General", "Delta_Time"}, 1.0, {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \optional_key{derivatives_mode_,Derivatives_Mode,string,"Covariant
   * Gaussian"}
   *
   * The mode of calculating the gradients, for example gradients of baryon
   * current. Currently SMASH supports two derivatives modes:
   *  - <tt>"Covariant Gaussian"</tt> and
   *  - <tt>"Finite difference"</tt>.
   *
   * Covariant Gaussian derivatives can be used when Covariant Gaussian smearing
   * is used; they are Lorentz covariant, but they do not calculate the time
   * derivative of the current properly. The `"Finite difference"` mode requires
   * using the lattice, and the derivatives are calculated based on finite
   * differences of a given quantity at adjacent lattice nodes; this mode is
   * numerically more efficient.
   */
  /**
   * \see_key{derivatives_mode_}
   */
  inline static const Key<DerivativesMode> gen_derivativesMode{
      {"General", "Derivatives_Mode"},
      DerivativesMode::CovariantGaussian,
      {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \optional_key{ensembles_,Ensembles,int,1}
   *
   * Number of parallel ensembles in the simulation.
   *
   * An ensemble is an instance of the system, and without mean-field potentials
   * it is practically equivalent to a completely separate and uncorrelated
   * event. Each ensemble is an independent simulation: initialization,
   * collisions, decays, box wall crossings, and propagation of particles is
   * performed independently within each ensemble.
   *
   * However, the densities and mean-field potentials are computed as averages
   * over all ensembles (within a given event). This process can be also viewed
   * as calculating densities and mean-fields by summing over particles in all
   * ensembles combined, where each particle contributes to the local charge
   * with a weight of 1/n_ensembles. Such technique is called the *parallel
   * ensemble* technique. It increases the statistics necessary for a precise
   * density calculation without increasing the number of collisions, which is
   * not the case in the *full ensemble* method (see <tt>\ref testparticles_
   * "Testparticles"</tt> description). Because of this, the parallel ensembles
   * technique is computationally faster than the full ensemble technique.
   */
  /**
   * \see_key{ensembles_}
   */
  inline static const Key<int> gen_ensembles{
      {"General", "Ensembles"}, 1, {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \optional_key{expansion_rate_,Expansion_Rate,double,0.1}
   *
   * Corresponds to the speed of expansion of the universe in non-Minkowski
   * metrics if <tt>\ref metric_type_ "Metric_Type"</tt> is any other than
   * `"NoExpansion"`.
   *
   * It corresponds to \f$b_r/l_0\f$ if the metric type is `"MasslessFRW"` or
   * `"MassiveFRW"`, and to the parameter b in the exponential expansion where
   * \f$a(t) ~ e^{bt/2}\f$.
   */
  /**
   * \see_key{expansion_rate_}
   */
  inline static const Key<double> gen_expansionRate{
      {"General", "Expansion_Rate"}, 0.1, {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \optional_key{metric_type_,Metric_Type,string,"NoExpansion"}
   *
   * Select which kind of expansion the metric should have. This needs only be
   * specified for the sphere modus. Possible values:
   * - `"NoExpansion"` &rarr; Default SMASH run, with Minkowski metric
   * - `"MasslessFRW"` &rarr; FRW expansion going as \f$t^{1/2}\f$
   * - `"MassiveFRW"` &rarr; FRW expansion going as \f$t^{2/3}\f$
   * - `"Exponential"` &rarr; FRW expansion going as \f$e^{t/2}\f$
   */
  /**
   * \see_key{metric_type_}
   */
  inline static const Key<ExpansionMode> gen_metricType{
      {"General", "Metric_Type"}, ExpansionMode::NoExpansion, {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \optional_key{rfdd_mode_,Rest_Frame_Density_Derivatives_Mode,string,"Off"}
   *
   * The mode of calculating the gradients of currents, decides whether the rest
   * frame density derivatives are computed (these derivatives are needed for
   * the VDF potentials, but not for the Skyrme potentials).
   */
  /**
   * \see_key{rfdd_mode_}
   */
  inline static const Key<RestFrameDensityDerivativesMode>
      gen_restFrameDensityDerivativeMode{
          {"General", "Rest_Frame_Density_Derivatives_Mode"},
          RestFrameDensityDerivativesMode::Off,
          {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \optional_key{smearing_mode_,Smearing_Mode,string,"Covariant Gaussian"}
   *
   * The mode of smearing for density calculation.
   *
   * Smearing is necessary to ensure a smooth gradient calculation, and it can
   * be thought of as smoothing out charge density fluctuations due to the
   * finite number of test-particles used. In general, this is done by
   * distributing the contribution to charge density from a given particle
   * according to some prescription. For example, in Gaussian smearing the
   * charge density of a particle is given by a Gaussian with some chosen width,
   * centered at the position of the particle; the Gaussian is normalized such
   * that integrating over the entire space yields the charge of the particle.
   * In result, the particle's charge is "smeared" over the space around it.
   * Note that the case with no smearing is recovered when the charge
   * contribution from each particle is taken to be a Dirac delta function
   * centered at the position of the particle.
   *
   * Currently, SMASH supports three smearing modes:
   * -# <tt>"Covariant Gaussian"</tt>\n
   *    This smearing represents the charge density of a particle as a Gaussian
   *    centered at the position of a particle; the user can specify the width
   *    and the cutoff of the Gaussian (the employed Gaussians, in principle
   *    non-zero over the entire available space, are "cut off" at some distance
   *    r_cut from the particle to improve calculation time). This smearing is
   *    Lorentz covariant which results in correct density profiles of
   *    relativistic systems. The downside of the smearing is its long
   *    computation time, as well as the fact that when the density is added to
   *    lattice nodes, it is done so by Euler approximation (using the density
   *    value at the lattice node), which does not conserve the number of
   *    particles on the lattice.
   * -# <tt>"Triangular"</tt>\n
   *    This smearing requires lattice; it represents the charge density of a
   *    particle in a given space direction as a "triangle" peaking at the
   *    particle's position and linearly decreasing over a specified range. The
   *    user specifies the range of the smearing in units of lattice spacings.
   *    This smearing is relatively fast, and it does conserve the number of
   *    particles on the lattice (due to the fact that the Euler integration is
   *    exact for a linear function).
   * -# <tt>"Discrete"</tt>\n
   *    This smearing requires lattice; the easiest of all smearing modes, it
   *    adds a specified portion of the particle's charge density to a node
   *    closest to the particle's position, and distributes the remainder evenly
   *    among the 6 nearest neighbor nodes. The user specifies the weight given
   *    to the center node; for example, if this weight is 1/3, then each of the
   *    six nearest neighbor nodes gets 1/9 of the particle's charge. This
   *    smearing is extremely fast, but is also rather coarse and requires using
   *    a large number of test-particles to produce smooth gradients.
   */
  /**
   * \see_key{smearing_mode_}
   */
  inline static const Key<SmearingMode> gen_smearingMode{
      {"General", "Smearing_Mode"}, SmearingMode::CovariantGaussian, {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \optional_key_no_line{gaussian_sigma_,Gaussian_Sigma,double,1.0}
   *
   * Parameter for Covariant Gaussian smearing: Width of Gaussian distributions
   * that represent Wigner density of particles, in fm.
   */
  /**
   * \see_key{gaussian_sigma_}
   */
  inline static const Key<double> gen_smearingGaussianSigma{
      {"General", "Gaussian_Sigma"}, 1.0, {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \optional_key_no_line{gauss_cutoff_in_sigma_,Gauss_Cutoff_In_Sigma,double,4.0}
   *
   * Parameter for Covariant Gaussian smearing: Distance in sigma at which
   * gaussian is considered 0.
   */
  /**
   * \see_key{gauss_cutoff_in_sigma_}
   */
  inline static const Key<double> gen_smearingGaussCutoffInSigma{
      {"General", "Gauss_Cutoff_In_Sigma"}, 4.0, {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \optional_key_no_line{triangular_range_,Triangular_Range,double,2.0}
   *
   * Parameter for Triangular smearing: Half of the base of a symmetric triangle
   * that represents particle density, in units of lattice spacings.
   */
  /**
   * \see_key{triangular_range_}
   */
  inline static const Key<double> gen_smearingTriangularRange{
      {"General", "Triangular_Range"}, 2.0, {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \optional_key_no_line{discrete_weight_,Discrete_Weight,double,0.333333}
   *
   * Parameter for Discrete smearing: Weight given to particle density at the
   * the center node; cannot be smaller than 1./7. (the boundary case of 1./7.
   * results in an even distribution of particle's density over the center node
   * and 6 neighboring nodes).
   */
  /**
   * \see_key{discrete_weight_}
   */
  inline static const Key<double> gen_smearingDiscreteWeight{
      {"General", "Discrete_Weight"}, 0.333333, {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \optional_key{testparticles_,Testparticles,int,1}
   *
   * Number of test-particles per real particle in the simulation.
   *
   * The number of initial sampled particles is increased by this factor,
   * while all cross sections are decreased by this factor. In this
   * way the mean free path does not change. Larger number of testparticles
   * helps to reduce spurious effects of geometric collision criterion
   * (see \iref{Cheng:2001dz}). It also reduces correlations related
   * to collisions and decays (but not the ones related to mean fields),
   * therefore the larger the number of testparticles, the closer the results
   * of the simulations should be to the solution of the Boltzmann equation.
   * These advantages come at a cost of a larger computational time.
   *
   * Testparticles are a way to increase statistics necessary for
   * precise density calculation, which is why they are needed for mean-field
   * potentials. The technique of using testparticles for mean field
   * is called the *full ensemble* technique. The number of collisions (and
   * consequently the simulation time) scales as square of the number of
   * testparticles, and that is why a full ensemble is slower than a parallel
   * ensemble.
   */
  /**
   * \see_key{testparticles_}
   */
  inline static const Key<int> gen_testparticles{
      {"General", "Testparticles"}, 1, {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \optional_key{time_step_mode_,Time_Step_Mode,string,"Fixed"}
   *
   * The mode of time stepping. Possible values:
   * - `"None"` &rarr; `Delta_Time` is set to the `End_Time`. This cannot be
   * used with potentials.
   * - `"Fixed"`&rarr; Fixed-sized time steps at which collision-finding grid is
   *   created. More efficient for systems with many particles. The `Delta_Time`
   *   is provided by user.
   *
   * For `Delta_Time` explanation see \ref delta_time_ "here".
   */
  /**
   * \see_key{time_step_mode_}
   */
  inline static const Key<TimeStepMode> gen_timeStepMode{
      {"General", "Time_Step_Mode"}, TimeStepMode::Fixed, {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \optional_key{use_grid_,Use_Grid,bool,true}
   *
   * - `true` &rarr; A grid is used to reduce the combinatorics of interaction
   * lookup.
   * - `false` &rarr; No grid is used.
   */
  /**
   * \see_key{use_grid_}
   */
  inline static const Key<bool> gen_useGrid{
      {"General", "Use_Grid"}, true, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_default_,default,string,ALL}
   *
   * It determines the default logging level for all areas
   */
  /**
   * \see_key{log_default_}
   */
  inline static const Key<einhard::LogLevel> log_default{
      {"Logging", "default"}, einhard::ALL, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_main_,Main,string,$\{default\}}
   */
  /**
   * \see_key{log_main_}
   */
  inline static const Key<einhard::LogLevel> log_main{{"Logging", "Main"},
                                                      {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_experiment_,Experiment,string,$\{default\}}
   */
  /**
   * \see_key{log_experiment_}
   */
  inline static const Key<einhard::LogLevel> log_experiment{
      {"Logging", "Experiment"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_box_,Box,string,$\{default\}}
   */
  /**
   * \see_key{log_box_}
   */
  inline static const Key<einhard::LogLevel> log_box{{"Logging", "Box"},
                                                     {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_collider_,Collider,string,$\{default\}}
   */
  /**
   * \see_key{log_collider_}
   */
  inline static const Key<einhard::LogLevel> log_collider{
      {"Logging", "Collider"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_sphere_,Sphere,string,$\{default\}}
   */
  /**
   * \see_key{log_sphere_}
   */
  inline static const Key<einhard::LogLevel> log_sphere{{"Logging", "Sphere"},
                                                        {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_action_,Action,string,$\{default\}}
   */
  /**
   * \see_key{log_action_}
   */
  inline static const Key<einhard::LogLevel> log_action{{"Logging", "Action"},
                                                        {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_input_parser_,InputParser,string,$\{default\}}
   */
  /**
   * \see_key{log_input_parser_}
   */
  inline static const Key<einhard::LogLevel> log_inputParser{
      {"Logging", "InputParser"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_particle_type_,ParticleType,string,$\{default\}}
   */
  /**
   * \see_key{log_particle_type_}
   */
  inline static const Key<einhard::LogLevel> log_particleType{
      {"Logging", "ParticleType"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_find_scatter_,FindScatter,string,$\{default\}}
   */
  /**
   * \see_key{log_find_scatter_}
   */
  inline static const Key<einhard::LogLevel> log_findScatter{
      {"Logging", "FindScatter"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_clock_,Clock,string,$\{default\}}
   */
  /**
   * \see_key{log_clock_}
   */
  inline static const Key<einhard::LogLevel> log_clock{{"Logging", "Clock"},
                                                       {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_decay_modes_,DecayModes,string,$\{default\}}
   */
  /**
   * \see_key{log_decay_modes_}
   */
  inline static const Key<einhard::LogLevel> log_decayModes{
      {"Logging", "DecayModes"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_resonances_,Resonances,string,$\{default\}}
   */
  /**
   * \see_key{log_resonances_}
   */
  inline static const Key<einhard::LogLevel> log_resonances{
      {"Logging", "Resonances"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_scatter_action_,ScatterAction,string,$\{default\}}
   */
  /**
   * \see_key{log_scatter_action_}
   */
  inline static const Key<einhard::LogLevel> log_scatterAction{
      {"Logging", "ScatterAction"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_distributions_,Distributions,string,$\{default\}}
   */
  /**
   * \see_key{log_distributions_}
   */
  inline static const Key<einhard::LogLevel> log_distributions{
      {"Logging", "Distributions"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_propagation_,Propagation,string,$\{default\}}
   */
  /**
   * \see_key{log_propagation_}
   */
  inline static const Key<einhard::LogLevel> log_propagation{
      {"Logging", "Propagation"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_grid_,Grid,string,$\{default\}}
   */
  /**
   * \see_key{log_grid_}
   */
  inline static const Key<einhard::LogLevel> log_grid{{"Logging", "Grid"},
                                                      {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_list_,List,string,$\{default\}}
   */
  /**
   * \see_key{log_list_}
   */
  inline static const Key<einhard::LogLevel> log_list{{"Logging", "List"},
                                                      {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_nucleus_,Nucleus,string,$\{default\}}
   */
  /**
   * \see_key{log_nucleus_}
   */
  inline static const Key<einhard::LogLevel> log_nucleus{{"Logging", "Nucleus"},
                                                         {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_density_,Density,string,$\{default\}}
   */
  /**
   * \see_key{log_density_}
   */
  inline static const Key<einhard::LogLevel> log_density{{"Logging", "Density"},
                                                         {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_pauli_blocking_,PauliBlocking,string,$\{default\}}
   */
  /**
   * \see_key{log_pauli_blocking_}
   */
  inline static const Key<einhard::LogLevel> log_pauliBlocking{
      {"Logging", "PauliBlocking"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_tmn_,Tmn,string,$\{default\}}
   */
  /**
   * \see_key{log_tmn_}
   */
  inline static const Key<einhard::LogLevel> log_tmn{{"Logging", "Tmn"},
                                                     {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_fpe_,Fpe,string,$\{default\}}
   */
  /**
   * \see_key{log_fpe_}
   */
  inline static const Key<einhard::LogLevel> log_fpe{{"Logging", "Fpe"},
                                                     {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_lattice_,Lattice,string,$\{default\}}
   */
  /**
   * \see_key{log_lattice_}
   */
  inline static const Key<einhard::LogLevel> log_lattice{{"Logging", "Lattice"},
                                                         {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_pythia_,Pythia,string,$\{default\}}
   */
  /**
   * \see_key{log_pythia_}
   */
  inline static const Key<einhard::LogLevel> log_pythia{{"Logging", "Pythia"},
                                                        {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_grandcan_thermalizer_,GrandcanThermalizer,string,$\{default\}}
   */
  /**
   * \see_key{log_grandcan_thermalizer_}
   */
  inline static const Key<einhard::LogLevel> log_grandcanThermalizer{
      {"Logging", "GrandcanThermalizer"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_cross_sections_,CrossSections,string,$\{default\}}
   */
  /**
   * \see_key{log_cross_sections_}
   */
  inline static const Key<einhard::LogLevel> log_crossSections{
      {"Logging", "CrossSections"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_output_,Output,string,$\{default\}}
   */
  /**
   * \see_key{log_output_}
   */
  inline static const Key<einhard::LogLevel> log_output{{"Logging", "Output"},
                                                        {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_hyper_surface_crossing_,HyperSurfaceCrossing,string,$\{default\}}
   */
  /**
   * \see_key{log_hyper_surface_crossing_}
   */
  inline static const Key<einhard::LogLevel> log_hyperSurfaceCrossing{
      {"Logging", "HyperSurfaceCrossing"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_initial_conditions_,InitialConditions,string,$\{default\}}
   */
  /**
   * \see_key{log_initial_conditions_}
   */
  inline static const Key<einhard::LogLevel> log_initialConditions{
      {"Logging", "InitialConditions"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_scatter_action_multi_,ScatterActionMulti,string,$\{default\}}
   */
  /**
   * \see_key{log_scatter_action_multi_}
   */
  inline static const Key<einhard::LogLevel> log_scatterActionMulti{
      {"Logging", "ScatterActionMulti"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{log_yaml_configuration_,YAML_Configuration,string,$\{default\}}
   */
  /**
   * \see_key{log_yaml_configuration_}
   */
  inline static const Key<einhard::LogLevel> log_yamlConfiguration{
      {"Logging", "YAML_Configuration"}, {"1.0"}};

  /// Alias for the type to be used in the list of keys.
  using key_references_variant = std::variant<
      std::reference_wrapper<const Key<bool>>,
      std::reference_wrapper<const Key<int>>,
      std::reference_wrapper<const Key<double>>,
      std::reference_wrapper<const Key<std::string>>,
      std::reference_wrapper<const Key<einhard::LogLevel>>,
      std::reference_wrapper<const Key<DerivativesMode>>,
      std::reference_wrapper<const Key<ExpansionMode>>,
      std::reference_wrapper<const Key<RestFrameDensityDerivativesMode>>,
      std::reference_wrapper<const Key<SmearingMode>>,
      std::reference_wrapper<const Key<TimeStepMode>>>;

  /// List of references to all existing SMASH keys.
  inline static const std::vector<key_references_variant> list = {
      std::cref(gen_endTime),
      std::cref(gen_modus),
      std::cref(gen_nevents),
      std::cref(gen_minNonEmptyEnsembles_number),
      std::cref(gen_minNonEmptyEnsembles_maximumEnsembles),
      std::cref(gen_randomseed),
      std::cref(gen_deltaTime),
      std::cref(gen_derivativesMode),
      std::cref(gen_ensembles),
      std::cref(gen_expansionRate),
      std::cref(gen_metricType),
      std::cref(gen_restFrameDensityDerivativeMode),
      std::cref(gen_smearingMode),
      std::cref(gen_smearingGaussianSigma),
      std::cref(gen_smearingGaussCutoffInSigma),
      std::cref(gen_smearingTriangularRange),
      std::cref(gen_smearingDiscreteWeight),
      std::cref(gen_testparticles),
      std::cref(gen_timeStepMode),
      std::cref(gen_useGrid),
      std::cref(log_default),
      std::cref(log_main),
      std::cref(log_experiment),
      std::cref(log_box),
      std::cref(log_collider),
      std::cref(log_sphere),
      std::cref(log_action),
      std::cref(log_inputParser),
      std::cref(log_particleType),
      std::cref(log_findScatter),
      std::cref(log_clock),
      std::cref(log_decayModes),
      std::cref(log_resonances),
      std::cref(log_scatterAction),
      std::cref(log_distributions),
      std::cref(log_propagation),
      std::cref(log_grid),
      std::cref(log_list),
      std::cref(log_nucleus),
      std::cref(log_density),
      std::cref(log_pauliBlocking),
      std::cref(log_tmn),
      std::cref(log_fpe),
      std::cref(log_lattice),
      std::cref(log_pythia),
      std::cref(log_grandcanThermalizer),
      std::cref(log_crossSections),
      std::cref(log_output),
      std::cref(log_hyperSurfaceCrossing),
      std::cref(log_initialConditions),
      std::cref(log_scatterActionMulti),
      std::cref(log_yamlConfiguration)};
};

}  // namespace smash

/*!\Userguide
 * \page minimum_nonempty_ensembles_
 * <hr>
 * \par Examples
 *
 * In the following example, the number of desired non-empty events is 1000
 * with a maximum number of 2000 events to be calculated. In this case the
 * calculation will stop either if 1000 events are not empty or 2000 events
 * have been calculated.
 * \verbatim
  General:
      Modus: Collider
      Minimum_Nonempty_Ensembles:
          Number: 1000
          Maximum_Ensembles_Run: 2000
      Ensembles: 1
  \endverbatim
 *
 * In contrast to the first example, in the next example we use 20 parallel
 * ensembles. Here, the maximum number of ensembles run is 2000. The calculation
 * will continue until either this number of ensembles is reached or 1000
 * ensembles contain interactions. Note that an event consists of 20 ensembles.
 * The 20 ensembles run in parallel, so the number of non-empty ensembles in the
 * ouput is between 1000 and 1019.
 * \verbatim
  General:
      Modus: Collider
      Minimum_Nonempty_Ensembles:
          Number: 1000
          Maximum_Ensembles_Run: 2000
      Ensembles: 20
  \endverbatim
 */

/*!\Userguide
 * \page input_logging_
 * <hr>
 * \par Example: Configuring the Logging Area
 *
 * To activate different logging levels for different logging areas, change the
 * default level for the desired areas. For example:
 *\verbatim
 Logging:
     default:    "WARN"
     Main:       "INFO"
     Experiment: "INFO"
     Pythia:     "DEBUG"
     Fpe:        "OFF"
 \endverbatim
 *
 * This will set all levels to `WARN` verbosity, still asking for informational
 * messages of `Main` and `Experiment` areas. Furthermore, `Pythia` debug
 * messages are requested, while any floating point exception message is turned
 * off.
 */

#endif  // SRC_INCLUDE_SMASH_VALIDATION_H_