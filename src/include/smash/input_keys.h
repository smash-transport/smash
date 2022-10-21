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
 * \li \subpage input_output_
 * \li \subpage input_lattice_
 * \li \subpage input_potentials_
 * \li \subpage input_forced_thermalization_
 *
 * \par Information on formatting of the input file
 *
 * The input file is made of sections, i.e. of keys containing as "value" a
 * series of keys and/or sections. In order to identify the content of a
 * section, it is important to keep a consistent indentation in the input file.
 * The convention is to use 4 spaces indentation in order to specify keys inside
 * a section.  For example:
 * \code
 * Output:
 *     Output_Interval: 1.0
 *     Particles:
 *         Format: ["Oscar2013"]
 * \endcode
 * This is a part of the input file. The `Output_Interval` key belongs to the
 * `Output` section, whereas `Particles` is in turn a section containing the
 * `Format` key.
 *
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
 * is generated. If the <tt>\ref key_gen_nevents_ "Nevents"</tt> key is not
 * specified, <b>this section with all its required keys must be present in the
 * SMASH input file</b>.
 *
 * number of ensembles is equal to the number of events, so that this option
 * will provide the desired number of non-empty events. \ref TBC
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

/*!\Userguide
 * \page input_collision_term_ Collision term
 *
 * The `Collision_Term` section in the input file can be used to configure SMASH
 * interactions. Before describing each possible key in detail, it is useful to
 * give some taste with a couple of examples.
 *
 * ### A real life example
 *
 * The following section in the input file configures SMASH to include all but
 * strangeness exchange involving 2 &harr; 2 scatterings, to treat N + Nbar
 * processes as resonance formations and to not force decays at the end of the
 * simulation. The elastic cross section is globally set to 30 mbarn and the
 * \f$ \sqrt{s} \f$ cutoff for elastic nucleon + nucleon collisions is 1.93 GeV.
 * All collisions are performed isotropically and 2 &harr; 1 processes are
 * forbidden.
 *
 *\verbatim
 Collision_Term:
     Included_2to2: ["Elastic","NN_to_NR","NN_to_DR","KN_to_KN","KN_to_KDelta"]
     Two_to_One: true
     Force_Decays_At_End: false
     NNbar_Treatment: "resonances"
     Elastic_Cross_Section: 30.0
     Elastic_NN_Cutoff_Sqrts: 1.93
     Isotropic: true
 \endverbatim
 *
 * If necessary, all collisions can be turned off by adding
 *\verbatim
     No_Collisions: True
 \endverbatim
 * in the configuration file.
 *
 * ### Configuring deuteron multi-particle reactions
 *
 * The following example configures SMASH to include deuteron multi-particle
 * reactions scatterings.
 *\verbatim
 Collision_Term:
     Collision_Criterion: Stochastic
     Multi_Particle_Reactions: ["Deuteron_3to2"]
 \endverbatim
 * Note, that the d' should not be included in the \e particels.txt file,
 * otherwise `PiDeuteron_to_pidprime` and `NDeuteron_to_Ndprime` have to be
 * excluded from `Included_2to2` by listing all 2-to-2 reactions except those
 * two.
 *
 * <hr>
 * In this page many generic keys are described. For information about further
 * tuning possibilities, see the following pages:
 * - \subpage input_collision_term_pauliblocker_
 * - \subpage input_collision_term_string_parameters_
 * - \subpage input_collision_term_dileptons_
 * - \subpage input_collision_term_photons_
 */

/*!\Userguide
 * \page input_collision_term_pauliblocker_ Pauli Blocking
 *
 * Pauli blocking can be activated and customized using the `Pauli_Blocking`
 * section within `Collision_Term`. For example:
 *\verbatim
 Collision_Term:
     Pauli_Blocking:
         Spatial_Averaging_Radius: 1.86
         Gaussian_Cutoff: 2.2
         Momentum_Averaging_Radius: 0.08
 \endverbatim
 */

/*!\Userguide
 * \page input_collision_term_string_parameters_ String Parameters
 *
 * Within `Collision_Term` section, the `String_Parameters` section can be used
 * to modify a series of parameters which affect the string fragmentation.
 */

/*!\Userguide
 * \page input_collision_term_dileptons_ Dileptons
 *
 * Dilepton production can be enabled in the corresponding `Dileptons`
 * section in the `Collision_Term` one of the configuration file.
 * Remember to also activate the dilepton output in the output section.
 */

/*!\Userguide
 * \page input_collision_term_photons_ Photons
 *
 * Photon production can be enabled in the corresponding `Photon` section
 * in the `Collision_Term` one of the configuration file.
 * Remember to also activate the photon output in the output section.
 */

/*!\Userguide
 * \page input_modi_ Modi
 *
 * The `Modi` section is the place where the specified <tt>\ref key_gen_modus_
 * "Modus"</tt> shall be configured. For each possibility refer to the
 * corresponding documentation page:
 * - \subpage input_modi_collider_
 * - \subpage input_modi_sphere_
 * - \subpage input_modi_box_
 * - \subpage input_modi_list_
 * - \subpage input_modi_listbox_
 *
 * The `Modi` section has to contain a section named after the chosen modus and
 * in it the corresponding customization takes place.
 */

/*!\Userguide
 * \page input_modi_collider_ Collider
 *
 * The `Collider` modus can be customized using the options here below.
 * To further configure the projectile, target and the impact parameter, see
 * - \subpage input_modi_collider_projectile_and_target_ and
 * - \subpage input_modi_collider_impact_parameter_.
 *
 * \attention
 * The incident energy can be specified in different ways and one (and only one)
 * of these must be used. Alternatively, one can specify the individual beam
 * energies or momenta in the `Projectile` and `Target` sections (see \ref
 * input_modi_collider_projectile_and_target_ for details). In this case,
 * one must give either `E_Tot` or `E_Kin` or `P_Lab` for both `Projectile`
 * and `Target`.
 *
 * <hr>
 */

/*!\Userguide
 * \page input_modi_collider_projectile_and_target_ Projectile and Target
 *
 * Within the `Collider` section, two sections can be used for further
 * customizations:
 * - `Projectile` &rarr; Section for projectile nucleus. The projectile will
 *   start at \f$z<0\f$ and fly in positive \f$z\f$-direction, at \f$x\ge 0\f$.
 * - `Target` &rarr; Section for target nucleus. The target will start at
 *   \f$z>0\f$ and fly in negative \f$z\f$-direction, at \f$x \le 0\f$.
 *
 * <b>All keys described here below can be specified in the `Projectile` and/or
 * in the `Target` section.</b> Examples are given after the keys description.
 */

/*!\Userguide
 * \page input_modi_collider_impact_parameter_ Impact Parameter
 *
 * Within the `Collider` section, the `Impact` section can be used to specify
 * information about the impact parameter, defined as the distance <b>in fm</b>
 * of the two straight lines that the center of masses of the nuclei travel on.
 * The separation of the two colliding nuclei is by default along the x-axis.
 *
 * \warning
 * Note that there are no safeguards to prevent you from specifying negative
 * impact parameters. The value chosen here is simply the x-component of
 * \f$\vec{\mkern0mu b}\f$. The result will be that the projectile and target
 * will have switched position in x.
 */

/*!\Userguide
 * \page input_modi_sphere_ Sphere
 */

/*!\Userguide
 * \page input_modi_box_ Box
 */

/*!\Userguide
 * \page input_modi_list_ List
 * The `List` modus provides a modus for hydro afterburner calculations. It
 * takes files with a list of particles in \ref oscar2013_format
 * "Oscar 2013 format" as an input. These particles are treated as a starting
 * setup. Multiple events per file are supported. In the following, the input
 * keys are listed with a short description, an example is given and some
 * information about the input particle files is provided.
 */

/*!\Userguide
 * \page input_modi_listbox_ ListBox
 *
 * The `ListBox` modus provides the possibility to initialize a box with a given
 * set of particles. This modus uses all functionality from the `List` modus
 * itself. The only difference is that one has to specify the length of the box.
 * Apart from that, the usage should be equivalent to \ref input_modi_list_
 * "the \c List modus". Refer to it for more details.
 *
 * ### Configuration example
 * \verbatim
 Modi:
     ListBox:
         File_Directory: "particle_lists_in"
         File_Prefix: "event"
         Shift_Id: 0
         Length: 10.0

 \endverbatim
 */

/*!\Userguide
 * \page input_output_ Output
 *
 * To produce a certain output content it is necessary to explicitly configure
 * it in the `Output` section of the configuration file. This means, that the
 * `Output` section needs to contain one or more subsection for each desired
 * content. Additionally, there are general output configuration parameters that
 * can be used for further customization.
 */

/*!\Userguide
 * \page input_lattice_ Lattice
 *
 * It is possible to configure a Lattice for the 3D space, which can be useful
 * to speed up the computation of the potentials. Note though, that this goes in
 * hand with a loss of accuracy: If the lattice is applied, the evaluation of
 * the potentials is carried out only on the nodes of the lattice. Intermediate
 * values are interpolated.
 *
 * The configuration of a lattice is usually not necessary, it is however
 * required if the \ref output_vtk_lattice_ "Thermodynamic VTK Output", the
 * \ref thermodyn_lattice_output_ "Thermodynamic Lattice Output" or the
 * <tt>\ref key_lattice_pot_affect_threshold_
 * "Potentials_Affect_Thresholds"</tt> option is enabled. To configure the
 * thermodynamic output, use \ref input_output_ "the \c Output section".
 *
 * The following parameters are only required, if the `Lattice` section is used
 * in the configuration file. Otherwise, no lattice will be used at all.
 */

/*!\Userguide
 * \page input_potentials_ %Potentials
 *
 * SMASH simulation supports two sets of potentials:
 * -# Skyrme and/or Symmetry potentials
 * -# VDF (vector density functional) model potentials,
 *    https://arxiv.org/pdf/2011.06635.pdf \n
 *
 * Coulomb potentials can be additionally enabled.
 *
 * Skyrme and VDF potentials both describe the behavior of symmetric nuclear
 * matter. The symmetry potential can adjust the Skyrme potential (but not the
 * VDF potential) to include effects due to isospin. The Skyrme and Symmetry
 * potentials are semi-relativistic, while the VDF potential is fully
 * relativistic.
 * - \subpage input_potentials_skyrme_
 * - \subpage input_potentials_symmetry_
 * - \subpage input_potentials_VDF_
 * - \subpage input_potentials_coulomb_
 *
 * ### Configuring Skyrme potentials
 *
 * The following snippet of the configuration file configures SMASH such
 * that the Skyrme as well as the Symmetry potential are activated for the
 * simulation. There is however no requirement to include both simultaneously.
 * They can be switched on and off individually.
 *\verbatim
 Potentials:
     Skyrme:
         Skyrme_A: -209.2
         Skyrme_B: 156.4
         Skyrme_Tau: 1.35
     Symmetry:
         S_Pot: 18.0
     Coulomb:
         R_Cut: 5.0
 \endverbatim
 * Note that the Coulomb potential requires a <tt>\ref input_lattice_
 * "Lattice"</tt> while for the other potentials it can be used as an
 * optimisation.
 *
 * ### Configuring VDF Potentials
 *
 * The following snippets from the configuration file configure SMASH such
 * that the VDF potential is activated for the simulation.
 *
 * In the first example, VDF potentials are configured to reproduce the default
 * SMASH Skyrme potentials (without the symmetry potential, as it is not
 * described within the VDF model):
 *\verbatim
 Potentials:
     VDF:
         Sat_rhoB: 0.168
         Powers: [2.0, 2.35]
         Coeffs: [-209.2, 156.5]
 \endverbatim
 *
 * In the second example, VDF potentials are configured to describe nuclear
 * matter with saturation density of \f$\rho_0 = \mathrm{0.160 fm}^{-3}\f$,
 * binding energy of \f$B_0 = -16.3\f$ MeV, the critical point of the
 * ordinary nuclear liquid-gas phase transition at \f$T_c^{(N)} = 18\f$ MeV and
 * \f$\rho_c^{(N)} = 0.375 \rho_0\f$, the critical point of the conjectured
 * "QGP-like" phase transition at \f$T_c^{(Q)} = 100\f$ MeV and
 * \f$\rho_c^{(Q)} = 3.0\rho_0\f$, and the boundaries of the spinodal region
 * of the "QGP-like" phase transition at \f$\eta_L = 2.50 \rho_0\f$ and
 * \f$\eta_R = 3.315 \rho_0\f$:
 *\verbatim
 Potentials:
     VDF:
         Sat_rhoB: 0.160
         Powers: [1.7681391, 3.5293515, 5.4352788, 6.3809822]
         Coeffs: [-8.450948e+01, 3.843139e+01, -7.958557e+00, 1.552594e+00]
 \endverbatim
 */

/*!\Userguide
 * \page input_potentials_skyrme_ Skyrme
 *
 * The Skyrme potential has the form
 * \f[ U_{Sk} = A(\rho/\rho_0) + B (\rho/\rho_0)^{\tau} \,, \f]
 * where \f$\rho\f$ is baryon density in the local Eckart rest frame.
 * Its parameters must be specified in the `Skyrme` subsection of the
 * `%Potentials` one.
 */

/*!\Userguide
 * \page input_potentials_symmetry_ Symmetry
 *
 * The symmetry potential has the form
 * \f[ U_{Sym} = \pm 2 S_{pot} \frac{I_3}{I} \frac{\rho_{I_3}}{\rho_0}
 * + S(\rho_B)\left(\frac{\rho_{I_3}}{\rho_B}\right)^2 \,, \f]
 * where \f$ \rho_{I_3}\f$ is the density of the relative isospin \f$ I_3/I
 * \f$ and \f$ \rho_B \f$ is the net baryon density and
 * \f[ S(\rho_B)=12.3\,\mathrm{MeV}\times
 * \left(\frac{\rho_B}{\rho_0}\right)^{2/3}+
 * 20\,\mathrm{MeV}\times\left(\frac{\rho_B}{\rho_0}\right)^\gamma\;. \f]
 * Parameters must be specified in the `Symmetry` subsection of the
 * `%Potentials` one.
 */

/*!\Userguide
 * \page input_potentials_VDF_ VDF
 *
 * The VDF potential is a four-vector of the form
 * \f[
 * A^{\mu} = \sum_{i=1}^N C_i
 * \left(\frac{\rho}{\rho_0}\right)^{b_i - 2}
 * \frac{j^{\mu}}{\rho_0} \,,
 * \f]
 * where \f$j^{\mu}\f$ is baryon 4-current, \f$\rho\f$ is baryon density in the
 * local Eckart rest frame, and \f$\rho_0\f$ is the saturation density. The
 * parameters of the potential, the coefficients \f$C_i\f$ and the powers
 * \f$b_i\f$, are fitted to reproduce a chosen set of properties of dense
 * nuclear matter, and in particular these may include describing two first
 * order phase transitions: the well-known phase transition in ordinary nuclear
 * matter, and a transition at high baryon densities meant to model a possible
 * QCD phase transition (a "QGP-like" phase transition); see
 * https://arxiv.org/pdf/2011.06635.pdf for details and example parameter sets
 * for the case \f$N=4\f$. The user can decide how many terms \f$N\f$ should
 * enter the potential by populating the coefficients and powers vectors in the
 * config file with a chosen number of entries. The number of coefficients must
 * match the number of powers.
 *
 * The potential parameters must be specified in the `VDF` subsection of the
 * `%Potentials` one.
 */

/*!\Userguide
 * \page input_potentials_coulomb_ Coulomb
 *
 * The Coulomb potential in SMASH includes the electric and magnetic field.
 * For simplicity we assume magnetostatics such that the fields can be
 * directly calculated as
 * \f[
 * \vec{E}(\vec{r})=-\vec{\nabla} \phi(\vec{r}) =
 * -\vec{\nabla}\int\frac{\rho(\vec{r}')}{|\vec{r}-\vec{r}'|} dV'
 * =\int\frac{\rho(\vec{r}')(\vec{r}-\vec{r}')}{|\vec{r}-\vec{r}'|^3}dV'
 * \f]
 * and
 * \f[
 * \vec{B}(\vec{r}) = \vec{\nabla}\times\vec{A}(\vec{r})
 * =\vec{\nabla}\times\int \frac{\vec{j}(\vec{r}')}{|\vec{r}-\vec{r}'|}dV'
 * =\int\vec{j}(\vec{r}')\times\frac{\vec{r}-\vec{r}'}{|\vec{r}-\vec{r}'|^3}dV'.
 * \f] These integrals are solved numerically on the SMASH lattice, where the
 * discretized equations read
 * \f[ \vec{E}(\vec{r}_j) = \sum_{i\neq j}
 * \frac{\rho(\vec{r}_i)(\vec{r}_j-\vec{r}_i)}{|\vec{r}_j-\vec{r}_i|^3}\Delta
 * V \f] and \f[ \vec{B}(\vec{r}_j)=\sum_{i\neq j} \vec{j}(\vec{r}_i)\times
 * \frac{\vec{r}_j-\vec{r}_i}{|\vec{r}_j-\vec{r}_i|^3} \Delta V \f]
 * with the lattice cell volume \f$ \Delta V \f$. For efficiency the integration
 * volume is cut at \f$ R_\mathrm{cut} \f$, which is taken from the
 * configuration. Note that in the final equations the summand for \f$i=j\f$
 * drops out because the contribution from that cell to the integral vanishes if
 * one assumes the current and density to be constant in the cell.
 */

/*!\Userguide
 * \page input_forced_thermalization_ Forced Thermalization
 *
 * Forced thermalization for certain regions is applied if the corresponding
 * `Forced_Thermalization` section is present in the configuration file.
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
   * ### Mandatory keys
   */

  /*!\Userguide
   * \page input_general_
   * \required_key_no_line{key_gen_end_time_,End_Time,double}
   *
   * The time <b>in fm</b> after which the evolution is stopped. Note
   * that the starting time depends on the chosen `Modus`.
   */
  /**
   * \see_key{key_gen_end_time_}
   */
  inline static const Key<double> gen_endTime{{"General", "End_Time"}, {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \required_key{key_gen_modus_,Modus,string}
   *
   * Selects a modus for the calculation, e.g.\ infinite matter
   * calculation, collision of two particles or collision of nuclei. The modus
   * will be configured in the <tt>\ref input_modi_ "Modi"</tt> section.
   * Recognized values are:
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
   * \see_key{key_gen_modus_}
   */
  inline static const Key<std::string> gen_modus{{"General", "Modus"}, {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \required_key{key_gen_nevents_,Nevents,int}
   *
   * Number of events to calculate.
   *
   * This key may be omitted on constraint that a minimum number
   * of ensembles containing interactions is requested, see
   * \subpage minimum_nonempty_ensembles_.
   */
  /**
   * \see_key{key_gen_nevents_}
   */
  inline static const Key<int> gen_nevents{{"General", "Nevents"}, {"1.0"}};

  /*!\Userguide
   * \page minimum_nonempty_ensembles_
   * \required_key{key_gen_mnee_number_,Number,int}
   *
   * The number of desired non-empty ensembles.\n
   */
  /**
   * \see_key{key_gen_mnee_number_}
   */
  inline static const Key<int> gen_minNonEmptyEnsembles_number{
      {"General", "Minimum_Nonempty_Ensembles", "Number"}, {"1.0"}};

  /*!\Userguide
   * \page minimum_nonempty_ensembles_
   * \required_key{key_gen_mnee_maximum_ensembles_,Maximum_Ensembles_Run,int}
   *
   * Maximum number of ensembles run. This number serves as a safeguard
   * against SMASH unexpectedly running for a long time.
   */
  /**
   * \see_key{key_gen_mnee_maximum_ensembles_}
   */
  inline static const Key<int> gen_minNonEmptyEnsembles_maximumEnsembles{
      {"General", "Minimum_Nonempty_Ensembles", "Maximum_Ensembles_Run"},
      {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \required_key{key_gen_randomseed_,Randomseed,int}
   *
   * Initial seed for the random number generator. If this is negative, the
   * seed will be randomly generated by the operating system.
   */
  /**
   * \see_key{key_gen_randomseed_}
   */
  inline static const Key<int> gen_randomseed{{"General", "Randomseed"},
                                              {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * <hr>
   * ### Optional keys
   */

  /*!\Userguide
   * \page input_general_
   * \optional_key_no_line{key_gen_delta_time_,Delta_Time,double,1.0}
   *
   * Fixed time step at which the collision-finding grid is recreated, and, if
   * potentials are on, momenta are updated according to the equations of
   * motion. The collision-finding grid finds all the collisions from time
   * t_{beginning_of_timestep} until time t_{beginning_of_timestep} +
   * `Delta_Time`, and puts them into a vector. The collisions are then sorted
   * in order of occurrence, and particles are propagated from collision to
   * collision. After each performed collision, additional collisions are found
   * for outgoing particles and merged into the sorted vector.
   *
   * If potentials are on, the Delta_Time should be small enough, typically
   * around 0.1 fm/c. However, if potentials are off, it can be arbitrarily
   * large. In this case it only influences the runtime, but not physics.
   * If `Time_Step_Mode = "None"` is chosen, then the user-provided value of
   * `Delta_Time` is ignored and `Delta_Time` is set to the `End_Time`.
   */
  /**
   * \see_key{key_gen_delta_time_}
   */
  inline static const Key<double> gen_deltaTime{
      {"General", "Delta_Time"}, 1.0, {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \optional_key{key_gen_derivatives_mode_,Derivatives_Mode,string,"Covariant
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
   * \see_key{key_gen_derivatives_mode_}
   */
  inline static const Key<DerivativesMode> gen_derivativesMode{
      {"General", "Derivatives_Mode"},
      DerivativesMode::CovariantGaussian,
      {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \optional_key{key_gen_ensembles_,Ensembles,int,1}
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
   * not the case in the *full ensemble* method (see <tt>\ref
   * key_gen_testparticles_ "Testparticles"</tt> description). Because of this,
   * the parallel ensembles technique is computationally faster than the full
   * ensemble technique.
   */
  /**
   * \see_key{key_gen_ensembles_}
   */
  inline static const Key<int> gen_ensembles{
      {"General", "Ensembles"}, 1, {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \optional_key{key_gen_expansion_rate_,Expansion_Rate,double,0.1}
   *
   * Corresponds to the speed of expansion of the universe in non-Minkowski
   * metrics if <tt>\ref key_gen_metric_type_ "Metric_Type"</tt> is any other
   * than `"NoExpansion"`.
   *
   * It corresponds to \f$b_r/l_0\f$ if the metric type is `"MasslessFRW"` or
   * `"MassiveFRW"`, and to the parameter b in the exponential expansion where
   * \f$a(t) ~ e^{bt/2}\f$.
   */
  /**
   * \see_key{key_gen_expansion_rate_}
   */
  inline static const Key<double> gen_expansionRate{
      {"General", "Expansion_Rate"}, 0.1, {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \optional_key{key_gen_metric_type_,Metric_Type,string,"NoExpansion"}
   *
   * Select which kind of expansion the metric should have. This needs only be
   * specified for the sphere modus. Possible values:
   * - `"NoExpansion"` &rarr; Default SMASH run, with Minkowski metric
   * - `"MasslessFRW"` &rarr; FRW expansion going as \f$t^{1/2}\f$
   * - `"MassiveFRW"` &rarr; FRW expansion going as \f$t^{2/3}\f$
   * - `"Exponential"` &rarr; FRW expansion going as \f$e^{t/2}\f$
   */
  /**
   * \see_key{key_gen_metric_type_}
   */
  inline static const Key<ExpansionMode> gen_metricType{
      {"General", "Metric_Type"}, ExpansionMode::NoExpansion, {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \optional_key{key_gen_rfdd_mode_,Rest_Frame_Density_Derivatives_Mode,string,"Off"}
   *
   * The mode of calculating the gradients of currents, decides whether the rest
   * frame density derivatives are computed (these derivatives are needed for
   * the VDF potentials, but not for the Skyrme potentials).
   */
  /**
   * \see_key{key_gen_rfdd_mode_}
   */
  inline static const Key<RestFrameDensityDerivativesMode>
      gen_restFrameDensityDerivativeMode{
          {"General", "Rest_Frame_Density_Derivatives_Mode"},
          RestFrameDensityDerivativesMode::Off,
          {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \optional_key{key_gen_smearing_mode_,Smearing_Mode,string,"Covariant
   * Gaussian"}
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
   * \see_key{key_gen_smearing_mode_}
   */
  inline static const Key<SmearingMode> gen_smearingMode{
      {"General", "Smearing_Mode"}, SmearingMode::CovariantGaussian, {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \optional_key_no_line{key_gen_gaussian_sigma_,Gaussian_Sigma,double,1.0}
   *
   * Parameter for Covariant Gaussian smearing: Width of Gaussian distributions
   * that represent Wigner density of particles, in fm.
   */
  /**
   * \see_key{key_gen_gaussian_sigma_}
   */
  inline static const Key<double> gen_smearingGaussianSigma{
      {"General", "Gaussian_Sigma"}, 1.0, {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \optional_key_no_line{key_gen_gauss_cutoff_in_sigma_,Gauss_Cutoff_In_Sigma,double,4.0}
   *
   * Parameter for Covariant Gaussian smearing: Distance in sigma at which
   * gaussian is considered 0.
   */
  /**
   * \see_key{key_gen_gauss_cutoff_in_sigma_}
   */
  inline static const Key<double> gen_smearingGaussCutoffInSigma{
      {"General", "Gauss_Cutoff_In_Sigma"}, 4.0, {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \optional_key_no_line{key_gen_triangular_range_,Triangular_Range,double,2.0}
   *
   * Parameter for Triangular smearing: Half of the base of a symmetric triangle
   * that represents particle density, in units of lattice spacings.
   */
  /**
   * \see_key{key_gen_triangular_range_}
   */
  inline static const Key<double> gen_smearingTriangularRange{
      {"General", "Triangular_Range"}, 2.0, {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \optional_key_no_line{key_gen_discrete_weight_,Discrete_Weight,double,0.333333}
   *
   * Parameter for Discrete smearing: Weight given to particle density at the
   * the center node; cannot be smaller than 1./7. (the boundary case of 1./7.
   * results in an even distribution of particle's density over the center node
   * and 6 neighboring nodes).
   */
  /**
   * \see_key{key_gen_discrete_weight_}
   */
  inline static const Key<double> gen_smearingDiscreteWeight{
      {"General", "Discrete_Weight"}, 0.333333, {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \optional_key{key_gen_testparticles_,Testparticles,int,1}
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
   * \see_key{key_gen_testparticles_}
   */
  inline static const Key<int> gen_testparticles{
      {"General", "Testparticles"}, 1, {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \optional_key{key_gen_time_step_mode_,Time_Step_Mode,string,"Fixed"}
   *
   * The mode of time stepping. Possible values:
   * - `"None"` &rarr; `Delta_Time` is set to the `End_Time`. This cannot be
   * used with potentials.
   * - `"Fixed"`&rarr; Fixed-sized time steps at which collision-finding grid is
   *   created. More efficient for systems with many particles. The `Delta_Time`
   *   is provided by user.
   *
   * For `Delta_Time` explanation see \ref key_gen_delta_time_ "here".
   */
  /**
   * \see_key{key_gen_time_step_mode_}
   */
  inline static const Key<TimeStepMode> gen_timeStepMode{
      {"General", "Time_Step_Mode"}, TimeStepMode::Fixed, {"1.0"}};

  /*!\Userguide
   * \page input_general_
   * \optional_key{key_gen_use_grid_,Use_Grid,bool,true}
   *
   * - `true` &rarr; A grid is used to reduce the combinatorics of interaction
   * lookup.
   * - `false` &rarr; No grid is used.
   */
  /**
   * \see_key{key_gen_use_grid_}
   */
  inline static const Key<bool> gen_useGrid{
      {"General", "Use_Grid"}, true, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_default_,default,string,ALL}
   *
   * It determines the default logging level for all areas
   */
  /**
   * \see_key{key_log_default_}
   */
  inline static const Key<einhard::LogLevel> log_default{
      {"Logging", "default"}, einhard::ALL, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_main_,Main,string,$\{default\}}
   */
  /**
   * \see_key{key_log_main_}
   */
  inline static const Key<einhard::LogLevel> log_main{{"Logging", "Main"},
                                                      {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_experiment_,Experiment,string,$\{default\}}
   */
  /**
   * \see_key{key_log_experiment_}
   */
  inline static const Key<einhard::LogLevel> log_experiment{
      {"Logging", "Experiment"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_box_,Box,string,$\{default\}}
   */
  /**
   * \see_key{key_log_box_}
   */
  inline static const Key<einhard::LogLevel> log_box{{"Logging", "Box"},
                                                     {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_collider_,Collider,string,$\{default\}}
   */
  /**
   * \see_key{key_log_collider_}
   */
  inline static const Key<einhard::LogLevel> log_collider{
      {"Logging", "Collider"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_sphere_,Sphere,string,$\{default\}}
   */
  /**
   * \see_key{key_log_sphere_}
   */
  inline static const Key<einhard::LogLevel> log_sphere{{"Logging", "Sphere"},
                                                        {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_action_,Action,string,$\{default\}}
   */
  /**
   * \see_key{key_log_action_}
   */
  inline static const Key<einhard::LogLevel> log_action{{"Logging", "Action"},
                                                        {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_input_parser_,InputParser,string,$\{default\}}
   */
  /**
   * \see_key{key_log_input_parser_}
   */
  inline static const Key<einhard::LogLevel> log_inputParser{
      {"Logging", "InputParser"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_particle_type_,ParticleType,string,$\{default\}}
   */
  /**
   * \see_key{key_log_particle_type_}
   */
  inline static const Key<einhard::LogLevel> log_particleType{
      {"Logging", "ParticleType"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_find_scatter_,FindScatter,string,$\{default\}}
   */
  /**
   * \see_key{key_log_find_scatter_}
   */
  inline static const Key<einhard::LogLevel> log_findScatter{
      {"Logging", "FindScatter"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_clock_,Clock,string,$\{default\}}
   */
  /**
   * \see_key{key_log_clock_}
   */
  inline static const Key<einhard::LogLevel> log_clock{{"Logging", "Clock"},
                                                       {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_decay_modes_,DecayModes,string,$\{default\}}
   */
  /**
   * \see_key{key_log_decay_modes_}
   */
  inline static const Key<einhard::LogLevel> log_decayModes{
      {"Logging", "DecayModes"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_resonances_,Resonances,string,$\{default\}}
   */
  /**
   * \see_key{key_log_resonances_}
   */
  inline static const Key<einhard::LogLevel> log_resonances{
      {"Logging", "Resonances"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_scatter_action_,ScatterAction,string,$\{default\}}
   */
  /**
   * \see_key{key_log_scatter_action_}
   */
  inline static const Key<einhard::LogLevel> log_scatterAction{
      {"Logging", "ScatterAction"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_distributions_,Distributions,string,$\{default\}}
   */
  /**
   * \see_key{key_log_distributions_}
   */
  inline static const Key<einhard::LogLevel> log_distributions{
      {"Logging", "Distributions"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_propagation_,Propagation,string,$\{default\}}
   */
  /**
   * \see_key{key_log_propagation_}
   */
  inline static const Key<einhard::LogLevel> log_propagation{
      {"Logging", "Propagation"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_grid_,Grid,string,$\{default\}}
   */
  /**
   * \see_key{key_log_grid_}
   */
  inline static const Key<einhard::LogLevel> log_grid{{"Logging", "Grid"},
                                                      {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_list_,List,string,$\{default\}}
   */
  /**
   * \see_key{key_log_list_}
   */
  inline static const Key<einhard::LogLevel> log_list{{"Logging", "List"},
                                                      {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_nucleus_,Nucleus,string,$\{default\}}
   */
  /**
   * \see_key{key_log_nucleus_}
   */
  inline static const Key<einhard::LogLevel> log_nucleus{{"Logging", "Nucleus"},
                                                         {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_density_,Density,string,$\{default\}}
   */
  /**
   * \see_key{key_log_density_}
   */
  inline static const Key<einhard::LogLevel> log_density{{"Logging", "Density"},
                                                         {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_pauli_blocking_,PauliBlocking,string,$\{default\}}
   */
  /**
   * \see_key{key_log_pauli_blocking_}
   */
  inline static const Key<einhard::LogLevel> log_pauliBlocking{
      {"Logging", "PauliBlocking"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_tmn_,Tmn,string,$\{default\}}
   */
  /**
   * \see_key{key_log_tmn_}
   */
  inline static const Key<einhard::LogLevel> log_tmn{{"Logging", "Tmn"},
                                                     {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_fpe_,Fpe,string,$\{default\}}
   */
  /**
   * \see_key{key_log_fpe_}
   */
  inline static const Key<einhard::LogLevel> log_fpe{{"Logging", "Fpe"},
                                                     {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_lattice_,Lattice,string,$\{default\}}
   */
  /**
   * \see_key{key_log_lattice_}
   */
  inline static const Key<einhard::LogLevel> log_lattice{{"Logging", "Lattice"},
                                                         {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_pythia_,Pythia,string,$\{default\}}
   */
  /**
   * \see_key{key_log_pythia_}
   */
  inline static const Key<einhard::LogLevel> log_pythia{{"Logging", "Pythia"},
                                                        {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_grandcan_thermalizer_,GrandcanThermalizer,string,$\{default\}}
   */
  /**
   * \see_key{key_log_grandcan_thermalizer_}
   */
  inline static const Key<einhard::LogLevel> log_grandcanThermalizer{
      {"Logging", "GrandcanThermalizer"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_cross_sections_,CrossSections,string,$\{default\}}
   */
  /**
   * \see_key{key_log_cross_sections_}
   */
  inline static const Key<einhard::LogLevel> log_crossSections{
      {"Logging", "CrossSections"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_output_,Output,string,$\{default\}}
   */
  /**
   * \see_key{key_log_output_}
   */
  inline static const Key<einhard::LogLevel> log_output{{"Logging", "Output"},
                                                        {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_hyper_surface_crossing_,HyperSurfaceCrossing,string,$\{default\}}
   */
  /**
   * \see_key{key_log_hyper_surface_crossing_}
   */
  inline static const Key<einhard::LogLevel> log_hyperSurfaceCrossing{
      {"Logging", "HyperSurfaceCrossing"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_initial_conditions_,InitialConditions,string,$\{default\}}
   */
  /**
   * \see_key{key_log_initial_conditions_}
   */
  inline static const Key<einhard::LogLevel> log_initialConditions{
      {"Logging", "InitialConditions"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_scatter_action_multi_,ScatterActionMulti,string,$\{default\}}
   */
  /**
   * \see_key{key_log_scatter_action_multi_}
   */
  inline static const Key<einhard::LogLevel> log_scatterActionMulti{
      {"Logging", "ScatterActionMulti"}, {"1.0"}};

  /*!\Userguide
   * \page input_logging_
   * \optional_key{key_log_yaml_configuration_,YAML_Configuration,string,$\{default\}}
   */
  /**
   * \see_key{key_log_yaml_configuration_}
   */
  inline static const Key<einhard::LogLevel> log_yamlConfiguration{
      {"Logging", "YAML_Configuration"}, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_
   * \optional_key{key_CT_two_to_one_,Two_to_One,bool,true}
   *
   * Enable 2 &harr; 1 processes (resonance formation and decays).
   */
  /**
   * \see_key{key_CT_two_to_one_}
   */
  inline static const Key<bool> collTerm_twoToOne{
      {"Collision_Term", "Two_to_One"}, true, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_
   * \optional_key{key_CT_included_2to2_,Included_2to2,list of strings,["All"]}
   *
   * List that contains all possible 2 &harr; 2 process categories. Each process
   * of the listed category can be performed within the simulation. Possible
   * categories are:
   * - `"Elastic"` &rarr; elastic binary scatterings
   * - `"NN_to_NR"` &rarr; nucleon + nucleon &harr; nucleon + resonance
   * - `"NN_to_DR"` &rarr; nucleon + nucleon &harr; delta + resonance
   * - `"KN_to_KN"` &rarr; kaon + nucleon &harr; kaon + nucleon
   * - `"KN_to_KDelta"` &rarr; kaon + nucleon &harr; kaon + delta
   * - `"Strangeness_exchange"` &rarr; processes with strangeness exchange
   * - `"NNbar"` &rarr; annihilation processes, when NNbar_treatment is set to
   *   resonances; this is superseded if NNbar_treatment is set to anything else
   * - `"PiDeuteron_to_NN"` &rarr; deuteron + pion &harr; nucleon + nucleon and
   *   its CPT-conjugate
   * - `"PiDeuteron_to_pidprime"` &rarr; deuteron + pion &harr; d' + pion
   * - `"NDeuteron_to_Ndprime"` &rarr; deuteron + (anti-)nucleon &harr;
   *   d' + (anti-)nucleon, and their CPT-conjugates
   * - `"All"` &rarr; include all binary processes, no necessity to list each
   *   single category
   *
   * Detailed balance is preserved by these reaction switches: if a forward
   * reaction is off then the reverse is automatically off too.
   */
  /**
   * \see_key{key_CT_included_2to2_}
   */
  inline static const Key<ReactionsBitSet> collTerm_includedTwoToTwo{
      {"Collision_Term", "Included_2to2"},
      ReactionsBitSet{}.set(),  // All interactions => all bit set
      {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_
   * \optional_key{key_CT_mp_reactions_,Multi_Particle_Reactions,list of
   * strings,[]}
   *
   * List of reactions with more than 2 in- or outgoing particles that contains
   * all possible multi-particle process categories. Multi particle reactions
   * only work with the stochastic collision criterion. Possible categories are:
   * - `"Meson_3to1"` &rarr; Mesonic 3-to-1 reactions:
   *   <table>
   *   <tr>
   *   <td> \f$\strut\pi^0\pi^+\pi^-\leftrightarrow\omega\f$
   *   <td> \f$\strut\pi^0\pi^+\pi^-\leftrightarrow\phi\f$
   *   <td> \f$\strut\eta\pi^+\pi^-\leftrightarrow\eta'\f$
   *   <td> \f$\strut\eta\pi^0\pi^0\leftrightarrow\eta'\f$
   *   </table>
   *   Since detailed balance is enforced, the corresponding decays also have to
   *   be added in decaymodes.txt to enable the reactions.
   * - `"Deuteron_3to2"` &rarr; Deuteron 3-to-2 reactions:
   *   <table>
   *   <tr>
   *   <td> \f$\strut \pi pn\leftrightarrow\pi d\f$
   *   <td> \f$\strut Npn\leftrightarrow Nd\f$
   *   <td> \f$\strut \bar{N}pn\leftrightarrow\bar{N}d\f$
   *   </table>
   *   The deuteron has to be uncommented in particles.txt as well.
   *   Do not uncomment d' or make sure to exclude 2-body reactions involving
   *   the d' (i.e. no `"PiDeuteron_to_pidprime"` and `"NDeuteron_to_Ndprime"`
   *   in `Included_2to2`). Otherwise, the deuteron reactions are implicitly
   *   double-counted.
   * - `"A3_Nuclei_4to2"` &rarr; Create or destroy A = 3 nuclei (triton, He-3,
   *   hypertriton) by 4 &harr; 2 catalysis reactions such as
   *   \f$X NNN \leftrightarrow X t\f$, where \f$X\f$ can be a pion, nucleon, or
   *   antinucleon.
   * - `"NNbar_5to2"` &rarr; 5-to-2 back-reaction for NNbar annihilation:
   *   \f$\pi^0\pi^+\pi^-\pi^+\pi^- \rightarrow N\bar{N}\f$. Since detailed
   *   balance is enforced, `NNbar_Treatment` has to be set to "two to five" for
   *   this option.
   */
  /**
   * \see_key{key_CT_mp_reactions_}
   */
  inline static const Key<MultiParticleReactionsBitSet>
      collTerm_multiParticleReactions{
          {"Collision_Term", "Multi_Particle_Reactions"},
          MultiParticleReactionsBitSet{}.reset(),  // Empty list => no bit set
          {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_
   * \optional_key{key_CT_force_decays_at_end_,Force_Decays_At_End,bool,true}
   *
   * - `true` &rarr; Force all resonances to decay after last timestep.
   * - `false` &rarr; Don't force decays (final output can contain resonances).
   */
  /**
   * \see_key{key_CT_force_decays_at_end_}
   */
  inline static const Key<bool> collTerm_forceDecaysAtEnd{
      {"Collision_Term", "Force_Decays_At_End"}, true, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_
   * \optional_key{key_CT_no_collisions_,No_Collisions,bool,false}
   *
   * Disable all possible collisions, only allow decays to occur if not
   * forbidden by other options. Useful for running SMASH as a decay
   * afterburner, but not recommended in general, because it breaks the detailed
   * balance.
   */
  /**
   * \see_key{key_CT_no_collisions_}
   */
  inline static const Key<bool> collTerm_noCollisions{
      {"Collision_Term", "No_Collisions"}, false, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_
   * \optional_key{key_CT_nnbar_treatment_,NNbar_Treatment,string,"strings"}
   *
   * - `"no annihilation"` &rarr; No annihilation of NNbar is performed.
   * - `"resonances"` &rarr; Annihilation through
   *   \f$N\bar{N}\rightarrow\rho h_1(1170)\f$; combined with
   *   \f$\rho\rightarrow\pi\pi\f$ and \f$h_1(1170)\rightarrow\pi\rho\f$, which
   *   gives 5 pions on average. This option requires `"NNbar"` to be enabled in
   *   <tt>\ref key_CT_included_2to2_ "Included_2to2"</tt>.
   * - `"two to five"` &rarr; Direct Annhilation of NNbar to \f$5\pi\f$,
   *   matching the resonance treatment:
   *   \f$N\bar{N}\rightarrow\pi^0\pi^+\pi^-\pi^+\pi^-\f$.
   *   This option requires `"NNbar_5to2"` to be enabled in
   *   <tt>\ref key_CT_mp_reactions_ "Multi_Particle_Reactions"</tt>.
   * - `"strings"` &rarr; Annihilation through string fragmentation.
   */
  /**
   * \see_key{key_CT_nnbar_treatment_}
   */
  inline static const Key<std::string> collTerm_nnbarTreatment{
      {"Collision_Term", "NNbar_Treatment"}, "strings", {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_
   * \optional_key{key_CT_use_aqm_,Use_AQM,bool,true}
   *
   * Turn on AQM cross-sections for exotic combination of particles
   * (baryon-baryon cross-sections are scaled from proton-proton high energy
   * parametrization, for example). This includes both elastic and non-elastic
   * contributions; non-elastic contributions go through string fragmentation.
   * Turning off strings or elastic collisions while leaving this on will
   * result in the corresponding part of the AQM cross-sections to also be off.
   * Cross-sections parametrization are scaled according to
   * \f[
   * \frac{\sigma^{\mathrm{AQM}}_{\mathrm{process}}}
   * {\sigma^{\mathrm{AQM}}_\mathrm{ref\_process}}
   * \sigma^{\mathrm{param}}_\mathrm{ref\_process}
   * \f]
   * where \f$ \sigma^{\mathrm{AQM}}_x = 40 \left( \frac{2}{3}
   * \right)^{n_\mathrm{meson}} (1 - 0.4 x^s_1) (1 - 0.4 x^s_2) \f$, with
   * \f$n_\mathrm{meson}\f$ being the number of mesons in the process,
   * \f$x^s_{1,2}\f$ the fraction of strange quarks in the participant.
   * "process" is then a generic process and "ref_process" a reference process
   * such as PP for which solid parametrizations exist.
   * (\iref{Bass:1998ca})
   */
  /**
   * \see_key{key_CT_use_aqm_}
   */
  inline static const Key<bool> collTerm_useAQM{
      {"Collision_Term", "Use_AQM"}, true, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_
   * \optional_key{key_CT_res_lifetime_mod_,Resonance_Lifetime_Modifier,double,1.0}
   *
   * Multiplicative factor by which to scale the resonance lifetimes up or down.
   * This additionally has the effect of modifying the initial densities by
   * the same factor in the case of a box initialized with thermal
   * multiplicities (see `Use_Thermal_Multiplicities` in \ref XXX and \ref YYY).
   *
   * \warning This option is not fully physically consistent with some of the
   * other assumptions used in SMASH; notably, modifying this value **will**
   * break detailed balance in any gas which allows resonances to collide
   * inelastically, as this option breaks the relationship between the width and
   * lifetime of resonances. Note as well that in such gases, using a value of
   * 0.0 is known to make SMASH hang; it is recommended to use a small non-zero
   * value instead in these cases.
   */
  /**
   * \see_key{key_CT_res_lifetime_mod_}
   */
  inline static const Key<double> collTerm_resonanceLifetimeModifier{
      {"Collision_Term", "Resonance_Lifetime_Modifier"}, 1.0, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_
   * \optional_key{key_CT_string_with_prob_,Strings_with_Probability,bool,true}
   *
   * - `true` &rarr;
   *   String processes are triggered according to a probability increasing
   *   smoothly with the collisional energy from 0 to 1 in a certain energy
   *   window. At energies beyond that window, all the inelastic scatterings are
   *   via strings, while at the energies below that window, all the scatterings
   *   are via non-string processes. One should be careful that in this
   *   approach, the scatterings via resoances are also suppressed in the
   *   intermediate energy region, and vanishes at high energies, e.g.
   *   \f$p\pi\rightarrow\Delta\rightarrow\Sigma K\f$
   *   can't happen at a collisional energy beyond 2.2 GeV in this approach.
   *   Therefore, the cross sections of the scatterings to the certain final
   *   states, which might be crucial for the production of the rare species,
   *   will be reduced at the high energies.
   * - `false` &rarr;
   *   String processes always happen as long as the collisional energy exceeds
   *   the threshold value by 0.9 GeV, and the parametrized total cross section
   *   is larger than the sum of cross sections contributed by the non-string
   *   processes. The string cross section is thus obtained by taking the
   *   difference between them.
   */
  /**
   * \see_key{key_CT_string_with_prob_}
   */
  inline static const Key<bool> collTerm_stringsWithProbability{
      {"Collision_Term", "Strings_with_Probability"}, true, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_
   * \optional_key{key_CT_elastic_cross_section_,Elastic_Cross_Section,double,-1.0}
   *
   * If a non-negative value is given, it will override the parametrized
   * elastic cross sections (which are energy-dependent) with a constant value
   * (in mb). This constant elastic cross section is used for all collisions.
   */
  /**
   * \see_key{key_CT_elastic_cross_section_}
   */
  inline static const Key<double> collTerm_elasticCrossSection{
      {"Collision_Term", "Elastic_Cross_Section"}, -1.0, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_
   * \optional_key{key_CT_isotropic_,Isotropic,bool,false}
   *
   * Do all collisions isotropically.
   */
  /**
   * \see_key{key_CT_isotropic_}
   */
  inline static const Key<bool> collTerm_isotropic{
      {"Collision_Term", "Isotropic"}, false, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_
   * \optional_key{key_CT_max_cs_,Maximum_Cross_Section,double,
   * 200</tt> or <tt>2000}
   *
   * The maximal cross section that should be used when looking for collisions.
   * This means that all particle pairs, whose transverse distance is smaller or
   * equal to \f$\sqrt{\sigma_\mathrm{max}/\pi}\f$, will be checked for
   * collisions. <b>The default value is usually set to 200 (in mb)</b> and this
   * value occurs in the Delta peak of the \f$\pi+p\f$ cross section. Many SMASH
   * cross sections diverge close at the threshold; these divergent parts are
   * effectively cut off. If deuteron production via d' is considered, then the
   * default is increased to 2000 mb to function correctly (see
   * \iref{Oliinychenko:2018ugs}). The maximal cross section is scaled with
   * <tt>\ref key_CT_cs_scaling_ "Cross_Section_Scaling"</tt> factor.
   */
  /**
   * \see_key{key_CT_max_cs_}
   */
  inline static const Key<double> collTerm_maximumCrossSection{
      {"Collision_Term", "Maximum_Cross_Section"}, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_
   * \optional_key{key_CT_fixed_min_cell_length_,Fixed_Min_Cell_Length,double,2.5}
   *
   * The (minimal) length <b>in fm</b> used for the grid cells of the stochastic
   * criterion, only. Collisions are searched within grid cells only. Cell
   * lengths are scaled up so that grid contains all particles if fraction of a
   * cell length would remain at end of the grid.
   */
  /**
   * \see_key{key_CT_fixed_min_cell_length_}
   */
  inline static const Key<double> collTerm_fixedMinCellLength{
      {"Collision_Term", "Fixed_Min_Cell_Length"}, 2.5, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_
   * \optional_key{key_CT_cs_scaling_,Cross_Section_Scaling,double,1.0}
   *
   * Scale all cross sections by a global factor.
   * \warning Most cross sections are constrained by experimental data. Scaling
   * them will therefore lead to nonphysical results and is only meant for
   * explorative studies.
   */
  /**
   * \see_key{key_CT_cs_scaling_}
   */
  inline static const Key<double> collTerm_crossSectionScaling{
      {"Collision_Term", "Cross_Section_Scaling"}, 1.0, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_
   * \optional_key{key_CT_additional_el_cs_,Additional_Elastic_Cross_Section,double,0.0}
   *
   * Add an additional constant contribution <b>in mb</b> to the elastic cross
   * section.
   * \warning Most elastic cross sections are constrained by experimental data.
   * Adding an additional contribution to them will therefore lead to
   * nonphysical results and is only meant for explorative studies.
   */
  /**
   * \see_key{key_CT_additional_el_cs_}
   */
  inline static const Key<double> collTerm_additionalElasticCrossSection{
      {"Collision_Term", "Additional_Elastic_Cross_Section"}, 0.0, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_
   * \optional_key{key_CT_include_decays_end_,Include_Weak_And_EM_Decays_At_The_End,bool,false}
   *
   * Enable to also perform weak and electro-magnetic decays at the end of the
   * simulation. If enabled, all decays in the *decaymodes.txt* file are
   * considered at the end, even for hadrons usually considered stable (i.e.
   * with an on-shell width larger than the width_cutoff), for example
   * \f$\Sigma\f$, \f$\pi\f$ or \f$\eta\f$. Note that for isospin violating
   * decay modes all possible isospin combination have to be manually specified
   * in the *decaymodes.txt* file.
   */
  /**
   * \see_key{key_CT_include_decays_end_}
   */
  inline static const Key<bool> collTerm_includeDecaysAtTheEnd{
      {"Collision_Term", "Include_Weak_And_EM_Decays_At_The_End"},
      false,
      {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_
   * \optional_key{key_CT_elastic_nn_cutoff_sqrts_,Elastic_NN_Cutoff_Sqrts,double,1.98}
   *
   * The elastic collisions between two nucleons with \f$\sqrt{s}\f$ below
   * the specified value (<b>in GeV</b>) cannot happen.
   * - `Elastic_NN_Cutoff_Sqrts` < 1.88 &rarr;
   *   Below the threshold energy of the elastic collision, no effect.
   * - `Elastic_NN_Cutoff_Sqrts` > 2.02 &rarr;
   *   Beyond the threshold energy of the inelastic collision
   *   \f$NN\rightarrow NN\pi\f$, not suggested.
   */
  /**
   * \see_key{key_CT_elastic_nn_cutoff_sqrts_}
   */
  inline static const Key<double> collTerm_elasticNNCutoffSqrts{
      {"Collision_Term", "Elastic_NN_Cutoff_Sqrts"}, 1.98, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_
   * \optional_key{key_CT_strings_,Strings,bool,
   * (\ref key_gen_modus_ "Modus"!="Box")}
   *
   * - `true` &rarr; String excitation is enabled
   * - `false` &rarr; String excitation is disabled
   */
  /**
   * \see_key{key_CT_strings_}
   */
  inline static const Key<bool> collTerm_strings{{"Collision_Term", "Strings"},
                                                 {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_
   * \optional_key{key_CT_collision_criterion_,Collision_Criterion,string,"Covariant"}
   *
   * The following collision criterions can be used.
   *
   * - `"Geometric"` &rarr; <b>Geometric collision criterion</b>\n
   *   The geometric collision criterion calculates the two-particle impact
   *   parameter as the closest approach distance in the two-particle
   *   center-of-momentum frame by boosting to the respective frame. The
   *   collision time used for the ordering is calculated as the time of the
   *   closest approach in the computational frame. For further details, see
   *   \iref{Bass:1998ca}.
   *
   * - `"Stochastic"` &rarr; <b>Stochastic collision criterion</b>\n
   *   The stochastic collision criterion employs a probability to decide
   *   whether particles collide inside a given space-time cell. The probability
   *   is derived directly from the scattering rate given by the Boltzmann
   *   equation. The stochastic criterion is the only criterion that allows to
   *   treat multi-particle reactions. For more details, see
   *   \iref{Staudenmaier:2021lrg}.
   *   \note
   *   The stochastic criterion is only applicable within limits. For example,
   *   it might not lead to reasonable results for very dilute systems like pp
   *   collisions. Futhermore, the fixed time step mode is required. The
   *   assumption for the criterion is that only one reaction per particle per
   *   timestep occurs. Therefore, small enough timesteps (<tt>\ref
   *   key_gen_delta_time_ "Delta_Time"</tt>) have to be used. In doubt, test if
   *   the results change with smaller timesteps. Since the probability value is
   *   not by defintion limited to 1 in case of large timesteps, an error is
   *   thrown if it gets larger than 1.
   *
   * - `"Covariant"` &rarr; <b>Covariant collision criterion</b>\n
   *   The covariant collision criterion uses a covariant expression of the
   *   two-particle impact parameter in the two-particle center-of-momentum
   *   frame, which allows for its calculation in the computational frame
   *   without boosting. Furthermore, it calculates the collision times used for
   *   the collision ordering in the two-particle center-of-momentum frame.
   *   Further details are described in \iref{Hirano:2012yy}.
   */
  /**
   * \see_key{key_CT_collision_criterion_}
   */
  inline static const Key<std::string> collTerm_collisionCriterion{
      {"Collision_Term", "Strings"}, "Covariant", {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_
   * \optional_key{key_CT_warn_high_prob_,Only_Warn_For_High_Probability,bool,false}
   *
   * Only warn and not error for reaction probabilities higher than 1.
   * This switch is meant for very long production runs with the stochastic
   * criterion. It has no effect on the other criteria. If enabled, it is user
   * responsibility to make sure that the warning, that the probability has
   * slipped above 1, is printed very rarely.
   */
  /**
   * \see_key{key_CT_warn_high_prob_}
   */
  inline static const Key<bool> collTerm_onlyWarnForHighProbability{
      {"Collision_Term", "Only_Warn_For_High_Probability"}, false, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_pauliblocker_
   * \optional_key{key_CT_PB_spatial_averaging_radius_,Spatial_Averaging_Radius,double,1.86}
   *
   * Radius <b>in fm</b> of sphere for averaging in the coordinate space.
   */
  /**
   * \see_key{key_CT_PB_spatial_averaging_radius_}
   */
  inline static const Key<double> collTerm_pauliBlocking_spatialAveragingRadius{
      {"Collision_Term", "Pauli_Blocking", "Spatial_Averaging_Radius"},
      1.86,
      {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_pauliblocker_
   * \optional_key{key_CT_PB_gaussian_cutoff_,Gaussian_Cutoff,double,2.2}
   *
   * Radius <b>in fm</b> at which Gaussians used for smoothing are cut.
   */
  /**
   * \see_key{key_CT_PB_gaussian_cutoff_}
   */
  inline static const Key<double> collTerm_pauliBlocking_gaussianCutoff{
      {"Collision_Term", "Pauli_Blocking", "Gaussian_Cutoff"}, 2.2, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_pauliblocker_
   * \optional_key{key_CT_PB_momentum_av_radius_,Momentum_Averaging_Radius,double,0.08}
   *
   * Radius <b>in GeV/c</b> of sphere for averaging in the momentum space.
   */
  /**
   * \see_key{key_CT_PB_momentum_av_radius_}
   */
  inline static const Key<double>
      collTerm_pauliBlocking_momentumAveragingRadius{
          {"Collision_Term", "Pauli_Blocking", "Momentum_Averaging_Radius"},
          0.08,
          {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_string_parameters_
   * \optional_key{key_CT_SP_string_tension_,String_Tension,double,1.0}
   *
   * String tension \f$\kappa\f$ <b>in GeV/fm</b> connecting massless quarks in
   * Hamiltonian, \f[H=|p_1|+|p_2|+\kappa |x_1-x_2|\;.\f]
   * This parameter is only used to determine particles' formation times
   * according to the yo-yo formalism (in the soft string routine for now).
   */
  /**
   * \see_key{key_CT_SP_string_tension_}
   */
  inline static const Key<double> collTerm_stringParam_stringTension{
      {"Collision_Term", "String_Parameters", "String_Tension"}, 1.0, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_string_parameters_
   * \optional_key{key_CT_SP_gluon_beta_,Gluon_Beta,double,0.5}
   *
   * Parameter \f$\beta\f$ in parton distribution function for gluons,
   * \f[\mathrm{PDF}_g(x) \propto \frac{1}{x}(1-x)^{\beta+1}\;.\f]
   */
  /**
   * \see_key{key_CT_SP_gluon_beta_}
   */
  inline static const Key<double> collTerm_stringParam_gluonBeta{
      {"Collision_Term", "String_Parameters", "Gluon_Beta"}, 0.5, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_string_parameters_
   * \optional_key{key_CT_SP_gluon_pmin_,Gluon_Pmin,double,0.001}
   *
   * Smallest possible scale for gluon lightcone momentum <b>in GeV</b>.
   * This is divided by \f$\sqrt{s}\f$ to get the minimum fraction to be sampled
   * from PDF shown in <tt>\ref key_CT_SP_gluon_beta_ "Gluon_Beta"</tt>.
   */
  /**
   * \see_key{key_CT_SP_gluon_pmin_}
   */
  inline static const Key<double> collTerm_stringParam_gluonPMin{
      {"Collision_Term", "String_Parameters", "Gluon_Pmin"}, 0.001, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_string_parameters_
   * \optional_key{key_CT_SP_quark_alpha_,Quark_Alpha,double,2.0}
   *
   * Parameter \f$\alpha\f$ in parton distribution function for quarks,
   * \f[\mathrm{PDF}_q\propto x^{\alpha-1}(1-x)^{\beta-1}\;.\f]
   */
  /**
   * \see_key{key_CT_SP_quark_alpha_}
   */
  inline static const Key<double> collTerm_stringParam_quarkAlpha{
      {"Collision_Term", "String_Parameters", "Quark_Alpha"}, 2.0, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_string_parameters_
   * \optional_key{key_CT_SP_quark_beta_,Quark_Beta,double,7.0}
   *
   * Parameter \f$\beta\f$ in PDF for quarks shown in <tt>\ref
   * key_CT_SP_quark_alpha_ "Quark_Alpha"</tt>.
   */
  /**
   * \see_key{key_CT_SP_quark_beta_}
   */
  inline static const Key<double> collTerm_stringParam_quarkBeta{
      {"Collision_Term", "String_Parameters", "Quark_Beta"}, 7.0, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_string_parameters_
   * \optional_key{key_CT_SP_strange_supp_,Strange_Supp,double,0.16}
   *
   * Strangeness suppression factor \f$\lambda\f$,
   * \f[\lambda=
   * \frac{P(s\bar{s})}{P(u\bar{u})\vphantom{\bar{d}}}=
   * \frac{P(s\bar{s})}{P(d\bar{d})}\;.
   * \f]
   * Defines the probability to produce a \f$s\bar{s}\f$ pair relative to
   * producing a light \f$q\bar{q}\f$ pair.
   */
  /**
   * \see_key{key_CT_SP_strange_supp_}
   */
  inline static const Key<double> collTerm_stringParam_strangeSuppression{
      {"Collision_Term", "String_Parameters", "Strange_Supp"}, 0.16, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_string_parameters_
   * \optional_key{key_CT_SP_diquark_supp_,Diquark_Supp,double,0.036}
   *
   * Diquark suppression factor. Defines the probability to produce a diquark
   * antidiquark pair relative to producing a qurk antiquark pair.
   */
  /**
   * \see_key{key_CT_SP_diquark_supp_}
   */
  inline static const Key<double> collTerm_stringParam_diquarkSuppression{
      {"Collision_Term", "String_Parameters", "Diquark_Supp"}, 0.036, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_string_parameters_
   * \optional_key{key_CT_SP_sigma_perp_,Sigma_Perp,double,0.42}
   *
   * Parameter \f$\sigma_\perp\f$ <b>in GeV</b> in the distribution for
   * transverse momentum transfer between colliding hadrons \f$p_\perp\f$ and
   * string mass \f$M_X\f$,
   * \f[
   * \frac{d^3N}{dM^2_Xd^2\mathbf{p_\perp}}\propto
   * \frac{1}{M_X^2} \exp\left(-\frac{p_\perp^2}{\sigma_\perp^2}\right)\;.
   * \f]
   */
  /**
   * \see_key{key_CT_SP_sigma_perp_}
   */
  inline static const Key<double> collTerm_stringParam_sigmaPerp{
      {"Collision_Term", "String_Parameters", "Sigma_Perp"}, 0.42, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_string_parameters_
   * \optional_key{key_CT_SP_stringz_a_,StringZ_A,double,2.0}
   *
   * Parameter \f$a\f$ in Pythia fragmentation function \f$f(z)\f$,
   * \f[f(z) = \frac{1}{z} (1-z)^a \exp\left(-b\frac{m_T^2}{z}\right)\;.\f]
   */
  /**
   * \see_key{key_CT_SP_stringz_a_}
   */
  inline static const Key<double> collTerm_stringParam_stringZA{
      {"Collision_Term", "String_Parameters", "StringZ_A"}, 2.0, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_string_parameters_
   * \optional_key{key_CT_SP_stringz_b_,StringZ_B,double,0.55}
   *
   * Parameter \f$b\f$ <b>in 1/GeV²</b> in Pythia fragmentation function shown
   * in <tt>\ref key_CT_SP_stringz_a_ "StringZ_A"</tt>.
   */
  /**
   * \see_key{key_CT_SP_stringz_b_}
   */
  inline static const Key<double> collTerm_stringParam_stringZB{
      {"Collision_Term", "String_Parameters", ""}, 0.55, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_string_parameters_
   * \optional_key{key_CT_SP_separate_fragment_bar_,Separate_Fragment_Baryon,bool,true}
   *
   * Whether to use a separate fragmentation function for leading baryons in
   * non-diffractive string processes.
   */
  /**
   * \see_key{key_CT_SP_separate_fragment_bar_}
   */
  inline static const Key<bool> collTerm_stringParam_separateFragmentBaryon{
      {"Collision_Term", "String_Parameters", "Separate_Fragment_Baryon"},
      true,
      {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_string_parameters_
   * \optional_key{key_CT_SP_stringz_a_leading_,StringZ_A_Leading,double,0.2}
   *
   * Parameter \f$a\f$ in Lund fragmentation function used to sample the light
   * cone momentum fraction of leading baryons in non-diffractive string
   * processes.
   */
  /**
   * \see_key{key_CT_SP_stringz_a_leading_}
   */
  inline static const Key<double> collTerm_stringParam_stringZALeading{
      {"Collision_Term", "String_Parameters", "StringZ_A_Leading"},
      0.2,
      {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_string_parameters_
   * \optional_key{key_CT_SP_stringz_b_leading_,StringZ_B_Leading,double,2.0}
   *
   * Parameter \f$b\f$ <b>in 1/GeV²</b> in Lund fraghmentation function used to
   * sample the light cone momentum fraction of leading baryons in
   * non-diffractive string processes.
   */
  /**
   * \see_key{key_CT_SP_stringz_b_leading_}
   */
  inline static const Key<double> collTerm_stringParam_stringZBLeading{
      {"Collision_Term", "String_Parameters", "StringZ_B_Leading"},
      2.0,
      {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_string_parameters_
   * \optional_key{key_CT_SP_string_sigma_t_,String_Sigma_T,double,0.5}
   *
   * Standard deviation <b>in GeV</b> in Gaussian for transverse momentum
   * distributed to string fragments during fragmentation.
   */
  /**
   * \see_key{key_CT_SP_string_sigma_t_}
   */
  inline static const Key<double> collTerm_stringParam_stringSigmaT{
      {"Collision_Term", "String_Parameters", "String_Sigma_T"}, 0.5, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_string_parameters_
   * \optional_key{key_CT_SP_form_time_factor_,Form_Time_Factor,double,1.0}
   *
   * Factor to be multiplied with the formation time of string fragments from
   * the soft string routine.
   */
  /**
   * \see_key{key_CT_SP_form_time_factor_}
   */
  inline static const Key<double> collTerm_stringParam_formTimeFactor{
      {"Collision_Term", "String_Parameters", "Form_Time_Factor"},
      1.0,
      {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_string_parameters_
   * \optional_key{key_CT_SP_power_part_formation_,Power_Particle_Formation,double,±1}
   *
   * The default value of this parameter is `+1` if
   * \f$\sqrt{s}<200\,\mathrm{GeV}\f$ and `-1` otherwise. If positive, the power
   * with which the cross section scaling factor of string fragments grows in
   * time until it reaches 1. If negative, the scaling factor will be constant
   * and jump to 1 once the particle forms.
   */
  /**
   * \see_key{key_CT_SP_power_part_formation_}
   */
  inline static const Key<double> collTerm_stringParam_powerParticleFormation{
      {"Collision_Term", "String_Parameters", "Power_Particle_Formation"},
      {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_string_parameters_
   * \optional_key{key_CT_SP_formation_time_,Formation_Time,double,1.0}
   *
   * Parameter for formation time in string fragmentation, <b>in fm</b>.
   */
  /**
   * \see_key{key_CT_SP_formation_time_}
   */
  inline static const Key<double> collTerm_stringParam_formationTime{
      {"Collision_Term", "String_Parameters", "Formation_Time"}, 1.0, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_string_parameters_
   * \optional_key{key_CT_SP_m_dependent_formation_t_,Mass_Dependent_Formation_Times,bool,false}
   *
   * Whether the formation time of string fragments should depend on their mass.
   * If it is set to `true`, the formation time is calculated as
   * \f$\tau = \sqrt{2}\frac{m}{\kappa} \f$.
   */
  /**
   * \see_key{key_CT_SP_m_dependent_formation_t_}
   */
  inline static const Key<bool> collTerm_stringParam_mDependentFormationTimes{
      {"Collision_Term", "String_Parameters", "Mass_Dependent_Formation_Times"},
      false,
      {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_string_parameters_
   * \optional_key{key_CT_SP_probability_p_to_duu_,Prob_proton_to_d_uu,double,1./3}
   *
   * Probability of splitting an (anti)nucleon into the quark it has only once
   * and the diquark it contains twice in terms of flavour in the soft string
   * routine.
   */
  /**
   * \see_key{key_CT_SP_probability_p_to_duu_}
   */
  inline static const Key<double> collTerm_stringParam_probabilityPToDUU{
      {"Collision_Term", "String_Parameters", "Prob_proton_to_d_uu"},
      1.0 / 3,
      {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_string_parameters_
   * \optional_key{key_CT_SP_popcorn_rate_,Popcorn_Rate,double,0.15}
   *
   * Parameter StringFlav:popcornRate, which determines production rate of
   * popcorn mesons in string fragmentation. It is possible to produce a popcorn
   * meson from the diquark end of a string with certain probability (i.e.,
   * diquark to meson + diquark).
   */
  /**
   * \see_key{key_CT_SP_popcorn_rate_}
   */
  inline static const Key<double> collTerm_stringParam_popcornRate{
      {"Collision_Term", "String_Parameters", "Popcorn_Rate"}, 0.15, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_dileptons_
   * \optional_key{key_CT_dileptons_decays_,Decays,bool,false}
   *
   * Whether or not to enable dilepton production from hadron decays.
   * This includes direct decays as well as Dalitz decays. Dilepton decays
   * additionally have to be uncommented in the used *decaymodes.txt* file
   * (see also \ref input_collision_term_dileptons_note_ "this note").
   */
  /**
   * \see_key{key_CT_dileptons_decays_}
   */
  inline static const Key<bool> collTerm_dileptons_decays{
      {"Collision_Term", "Dileptons", "Decays"}, false, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_photons_
   * \required_key{key_CT_photons_fractional_photons,Fractional_Photons,int}
   *
   * Number of fractional photons sampled per single perturbatively produced
   * photon.
   */
  /**
   * \see_key{key_CT_photons_fractional_photons}
   */
  inline static const Key<int> collTerm_photons_fractionalPhotons{
      {"Collision_Term", "Photons", "Fractional_Photons"}, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_photons_
   * \optional_key{key_CT_photons_2to2_scatterings_,2to2_Scatterings,bool,false}
   *
   * Whether or not to enable photon production in mesonic scattering processes.
   */
  /**
   * \see_key{key_CT_photons_2to2_scatterings_}
   */
  inline static const Key<bool> collTerm_photons_twoToTwoScatterings{
      {"Collision_Term", "Photons", "2to2_Scatterings"}, false, {"1.0"}};

  /*!\Userguide
   * \page input_collision_term_photons_
   * \optional_key{key_CT_photons_bremsstrahlung_,Bremsstrahlung,bool,false}
   *
   * Whether or not to enable photon production in bremsstrahlung processes.
   */
  /**
   * \see_key{key_CT_photons_bremsstrahlung_}
   */
  inline static const Key<bool> collTerm_photons_bremsstrahlung{
      {"Collision_Term", "Photons", "Bremsstrahlung"}, false, {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_
   *
   * \par Ways to specify incident energies &rarr; Only one can be given!
   *
   * \required_key_no_line{key_MC_sqrtsnn_,Sqrtsnn,double}
   *
   * Defines the energy of the collision as center-of-mass energy in the
   * collision of one participant each from both nuclei (using the average
   * participant mass in the given nucleus). This key can be omitted if
   * the incident energy is specified in a different way.
   */
  /**
   * \see_key{key_MC_sqrtsnn_}
   */
  inline static const Key<double> modi_collider_sqrtSNN{
      {"Modi", "Collider", "Sqrtsnn"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_
   * \required_key_no_line{key_MC_e_kin_,E_Kin,double}
   *
   * Defines the energy of the collision by the kinetic energy per nucleon of
   * the projectile nucleus, <b>in AGeV</b>. This assumes the target nucleus is
   * at rest. Note, this can also be given per-beam as described in \ref
   * input_modi_collider_projectile_and_target_. This key can be omitted if
   * the incident energy is specified in a different way.
   */
  /**
   * \see_key{key_MC_e_kin_}
   */
  inline static const Key<double> modi_collider_eKin{
      {"Modi", "Collider", "E_Kin"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_
   * \required_key_no_line{key_MC_e_tot_,E_Tot,double}
   *
   * Defines the energy of the collision by the total energy per nucleon of
   * the projectile nucleus, <b>in AGeV</b>. This assumes the target nucleus is
   * at rest. Note, this can also be given per-beam as described in \ref
   * input_modi_collider_projectile_and_target_. This key can be omitted if
   * the incident energy is specified in a different way.
   */
  /**
   * \see_key{key_MC_e_tot_}
   */
  inline static const Key<double> modi_collider_eTot{
      {"Modi", "Collider", "E_Tot"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_
   * \required_key_no_line{key_MC_p_lab_,P_Lab,double}
   *
   * Defines the energy of the collision by the initial momentum per nucleon
   * of the projectile nucleus, <b>in AGeV</b>. This assumes the target nucleus
   * is at rest.  This must be positive.  Note, this can also be given per-beam
   * as described in \ref input_modi_collider_projectile_and_target_. This key
   * can be omitted if the incident energy is specified in a different way.
   */
  /**
   * \see_key{key_MC_p_lab_}
   */
  inline static const Key<double> modi_collider_pLab{
      {"Modi", "Collider", "P_Lab"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_
   * \optional_key{key_MC_calc_frame_,Calculation_Frame,string,"center of
   * velocity"}
   *
   * The frame in which the collision is calculated. Possible values are
   * - `"center of velocity"`
   * - `"center of mass"`
   * - `"fixed target"`
   *
   * \note
   * Using `E_Tot`, `E_kin` or `P_Lab` to quantify the collision energy is not
   * sufficient to configure a collision in a fixed target frame. You need to
   * additionally change the `Calculation_Frame`. Any format of incident energy
   * can however be combined with any calculation frame, the provided incident
   * energy is then intrinsically translated to the quantity needed for the
   * computation.
   */
  /**
   * \see_key{key_MC_calc_frame_}
   */
  inline static const Key<std::string> modi_collider_calculationFrame{
      {"Modi", "Collider", "Calculation_Frame"}, "center of velocity", {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_
   * \optional_key{key_MC_fermi_motion_,Fermi_Motion,string,"off"}
   *
   * - `"on"` &rarr; Switch Fermi motion on, it is recommended to also activate
   * potentials.
   * - `"off"` &rarr; Switch Fermi motion off.
   * - `"frozen"` &rarr; Use "frozen" if you want to use Fermi motion
   * without potentials.
   */
  /**
   * \see_key{key_MC_fermi_motion_}
   */
  inline static const Key<std::string> modi_collider_fermiMotion{
      {"Modi", "Collider", "Fermi_Motion"}, "off", {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_
   * \optional_key{key_MC_collision_within_nucleus_,Collisions_Within_Nucleus,bool,false}
   *
   * Determine whether to allow the first collisions within the same nucleus.
   * - `true` &rarr; First collisions within the same nucleus allowed.
   * - `false` &rarr; First collisions within the same nucleus forbidden.
   */
  /**
   * \see_key{key_MC_collision_within_nucleus_}
   */
  inline static const Key<bool> modi_collider_collisionWithinNucleus{
      {"Modi", "Collider", "Collisions_Within_Nucleus"}, false, {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_
   * \optional_key{key_MC_initial_distance_,Initial_Distance,double,2.0}
   *
   * The initial distance of the two nuclei <b>in fm</b>:
   * \f$z_{\rm min}^{\rm target} - z_{\rm max}^{\rm projectile}\f$.
   *
   * Note that this distance is applied before the Lorentz boost to the chosen
   * calculation frame, and thus the actual distance may be different.
   */
  /**
   * \see_key{key_MC_initial_distance_}
   */
  inline static const Key<double> modi_collider_initialDistance{
      {"Modi", "Collider", "Initial_Distance"}, 2.0, {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_projectile_and_target_
   * \required_key{key_MC_PT_particles_,%Particles,map<int\,int>}
   *
   * A map in which the keys are PDG codes and the values are number of
   * particles with that PDG code that should be in the current nucleus.
   * For example:
   * - `{2212: 82, 2112: 126}` &rarr; a lead-208 nucleus (82 protons and 126
   *   neutrons = 208 nucleons)
   * - `{2212: 1, 2112: 1, 3122: 1}` &rarr; for Hyper-Triton (one proton, one
   *   neutron and one \f$\Lambda\f$).
   */
  /**
   * \see_key{key_MC_PT_particles_}
   */
  inline static const Key<std::map<PdgCode, int>>
      modi_collider_projectile_particles{
          {"Modi", "Collider", "Projectile", "Particles"}, {"1.0"}};
  /**
   * \see_key{key_MC_PT_particles_}
   */
  inline static const Key<std::map<PdgCode, int>>
      modi_collider_target_particles{
          {"Modi", "Collider", "Target", "Particles"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_projectile_and_target_
   * \optional_key{key_MC_PT_diffusiveness_,Diffusiveness,double,</tt>\f$d(A)\f$<tt>}
   *
   * Diffusiveness of the Woods-Saxon distribution for the nucleus <b>in fm</b>.
   * In general, the default value is
   * \f[
   * d(A)=\begin{cases}
   * 0.545 & A \le 16\\
   * 0.54  & A > 16
   * \end{cases}\;.
   * \f]
   * For copper, zirconium, ruthenium, xenon, gold, lead and uranium, a more
   * specific default value is used (see nucleus.cc).
   */
  /**
   * \see_key{key_MC_PT_diffusiveness_}
   */
  inline static const Key<double> modi_collider_projectile_diffusiveness{
      {"Modi", "Collider", "Projectile", "Diffusiveness"}, {"1.0"}};
  /**
   * \see_key{key_MC_PT_diffusiveness_}
   */
  inline static const Key<double> modi_collider_target_diffusiveness{
      {"Modi", "Collider", "Target", "Diffusiveness"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_projectile_and_target_
   * \optional_key{key_MC_PT_radius_,Radius,double,</tt>\f$r(A)\f$<tt>}
   *
   * Radius of nucleus <b>in fm</b>. In general, the default value is
   * \f[
   * r(A)=\begin{cases}
   * 1.2  \, A^{1/3}                     & A \le 16\\
   * 1.12 \, A^{1/3} - 0.86 \, A^{-1/3}  & A > 16
   * \end{cases}\;.
   * \f]
   * For copper, zirconium, ruthenium, xenon, gold, lead, and uranium, a more
   * specific default value is used (see nucleus.cc).
   */
  /**
   * \see_key{key_MC_PT_radius_}
   */
  inline static const Key<double> modi_collider_projectile_radius{
      {"Modi", "Collider", "Projectile", "Radius"}, {"1.0"}};
  /**
   * \see_key{key_MC_PT_radius_}
   */
  inline static const Key<double> modi_collider_target_radius{
      {"Modi", "Collider", "Target", "Radius"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_projectile_and_target_
   * \optional_key{key_MC_PT_saturation_density_,Saturation_Density,double,
   * </tt>\f$\int\rho(r)\:\mathrm{d}^3r=N_{nucleons}\f$<tt>}
   *
   * Saturation density of the nucleus <b>in 1/fm³</b>.
   * If not any value is specified, the saturation density is calculated such
   * that the integral over the Woods-Saxon distribution returns the number of
   * nucleons in the nucleus.
   */
  /**
   * \see_key{key_MC_PT_saturation_density_}
   */
  inline static const Key<double> modi_collider_projectile_saturationDensity{
      {"Modi", "Collider", "Projectile", "Saturation_Density"}, {"1.0"}};
  /**
   * \see_key{key_MC_PT_saturation_density_}
   */
  inline static const Key<double> modi_collider_target_saturationDensity{
      {"Modi", "Collider", "Target", "Saturation_Density"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_projectile_and_target_
   * <hr>
   * \par Possible incident energies given per beam
   *
   * \required_key_no_line{key_MC_PT_e_tot_,E_Tot,double}
   *
   * Set the totat energy <b>in GeV</b> per particle of the beam. This key,
   * if used, must be present in both `Projectile` and `Target` section. This
   * key can be omitted if the incident energy is specified in a different way.
   */
  /**
   * \see_key{key_MC_PT_e_tot_}
   */
  inline static const Key<double> modi_collider_projectile_eTot{
      {"Modi", "Collider", "Projectile", "E_Tot"}, {"1.0"}};
  /**
   * \see_key{key_MC_PT_e_tot_}
   */
  inline static const Key<double> modi_collider_target_eTot{
      {"Modi", "Collider", "Target", "E_Tot"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_projectile_and_target_
   * \required_key_no_line{key_MC_PT_e_kin_,E_Kin,double}
   *
   * Set the kinetic energy <b>in GeV</b> per particle of the beam. This key,
   * if used, must be present in both `Projectile` and `Target` section. This
   * key can be omitted if the incident energy is specified in a different way.
   */
  /**
   * \see_key{key_MC_PT_e_kin_}
   */
  inline static const Key<double> modi_collider_projectile_eKin{
      {"Modi", "Collider", "Projectile", "E_Kin"}, {"1.0"}};
  /**
   * \see_key{key_MC_PT_e_kin_}
   */
  inline static const Key<double> modi_collider_target_eKin{
      {"Modi", "Collider", "Target", "E_Kin"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_projectile_and_target_
   * \required_key_no_line{key_MC_PT_p_lab_,P_Lab,double}
   *
   * Set the momentum <b>in GeV/c</b> per particle of the beam. This key,
   * if used, must be present in both `Projectile` and `Target` section. This
   * key can be omitted if the incident energy is specified in a different way.
   *
   * \note
   * If the beam specific kinetic energy or momentum is set using either of
   * these keys, then it must be specified in the same way (not necessarily same
   * value) for both beams. This is for example useful to simulate for p-Pb
   * collisions at the LHC, where the centre-of-mass system does not correspond
   * to the laboratory system (see \ref
   * input_modi_collider_projectile_and_target_ex1_ "example").
   */
  /**
   * \see_key{key_MC_PT_p_lab_}
   */
  inline static const Key<double> modi_collider_projectile_pLab{
      {"Modi", "Collider", "Projectile", "P_Lab"}, {"1.0"}};
  /**
   * \see_key{key_MC_PT_p_lab_}
   */
  inline static const Key<double> modi_collider_target_pLab{
      {"Modi", "Collider", "Target", "P_Lab"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_projectile_and_target_
   * <hr>
   * ### Custom nuclei
   *
   * It is possible to further customize the projectile and/or target using the
   * `Custom` section, which should then contain few required keys, if given.
   *
   * \required_key_no_line{key_MC_PT_custom_file_dir_,File_Directory,string}
   *
   * The directory where the external list with the nucleon configurations
   * is located. <b>Make sure to use an absolute path!</b>
   */
  /**
   * \see_key{key_MC_PT_custom_file_dir_}
   */
  inline static const Key<std::string>
      modi_collider_projectile_custom_fileDirectory{
          {"Modi", "Collider", "Projectile", "Custom", "File_Directory"},
          {"1.0"}};
  /**
   * \see_key{key_MC_PT_custom_file_dir_}
   */
  inline static const Key<std::string>
      modi_collider_target_custom_fileDirectory{
          {"Modi", "Collider", "Target", "Custom", "File_Directory"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_projectile_and_target_
   * \required_key_no_line{key_MC_PT_custom_file_name_,File_Name,string}
   *
   * The file name of the external list with the nucleon configurations.
   */
  /**
   * \see_key{key_MC_PT_custom_file_name_}
   */
  inline static const Key<std::string> modi_collider_projectile_custom_fileName{
      {"Modi", "Collider", "Projectile", "Custom", "File_Name"}, {"1.0"}};
  /**
   * \see_key{key_MC_PT_custom_file_name_}
   */
  inline static const Key<std::string> modi_collider_target_custom_fileName{
      {"Modi", "Collider", "Target", "Custom", "File_Name"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_projectile_and_target_
   * <hr>
   * ### Deformed nuclei
   *
   * It is possible to deform the projectile and/or target nuclei using the
   * `Deformed` section, which should then contain some configuration, if given.
   *
   * \required_key_no_line{key_MC_PT_deformed_auto_,Automatic,bool}
   *
   * - `true` &rarr; Set parameters of spherical deformation based on mass
   *   number of the nucleus. Currently the following deformed nuclei are
   *   implemented: Cu, Zr, Ru, Au, Pb, U, and Xe (see deformednucleus.cc). If
   *   set to true the other parameters should not be provided.
   * - `false` &rarr; Manually set parameters of spherical deformation. This
   *   requires the additional specification of `Beta_2`, `Beta_3`, `Beta_4`,
   *   which follow \iref{Moller:1993ed} and \iref{Schenke:2019ruo}.
   */
  /**
   * \see_key{key_MC_PT_deformed_auto_}
   */
  inline static const Key<bool> modi_collider_projectile_deformed_automatic{
      {"Modi", "Collider", "Projectile", "Deformed", "Automatic"}, {"1.0"}};
  /**
   * \see_key{key_MC_PT_deformed_auto_}
   */
  inline static const Key<bool> modi_collider_target_deformed_automatic{
      {"Modi", "Collider", "Target", "Deformed", "Automatic"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_projectile_and_target_
   * \required_key_no_line{key_MC_PT_deformed_betaII_,Beta_2,double}
   *
   * The deformation coefficient for the spherical harmonic Y_2_0 in the beta
   * decomposition of the nuclear radius in the deformed Woods-Saxon
   * distribution. This key can be omitted for automatic deformation
   * (`Automatic: true`).
   */
  /**
   * \see_key{key_MC_PT_deformed_betaII_}
   */
  inline static const Key<double> modi_collider_projectile_deformed_beta2{
      {"Modi", "Collider", "Projectile", "Deformed", "Beta_2"}, {"1.0"}};
  /**
   * \see_key{key_MC_PT_deformed_betaII_}
   */
  inline static const Key<double> modi_collider_target_deformed_beta2{
      {"Modi", "Collider", "Target", "Deformed", "Beta_2"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_projectile_and_target_
   * \required_key_no_line{key_MC_PT_deformed_betaIII_,Beta_3,double}
   *
   * The deformation coefficient for the spherical harmonic Y_3_0. This key can
   * be omitted for automatic deformation (`Automatic: true`).
   */
  /**
   * \see_key{key_MC_PT_deformed_betaIII_}
   */
  inline static const Key<double> modi_collider_projectile_deformed_beta3{
      {"Modi", "Collider", "Projectile", "Deformed", "Beta_3"}, {"1.0"}};
  /**
   * \see_key{key_MC_PT_deformed_betaIII_}
   */
  inline static const Key<double> modi_collider_target_deformed_beta3{
      {"Modi", "Collider", "Target", "Deformed", "Beta_3"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_projectile_and_target_
   * \required_key_no_line{key_MC_PT_deformed_betaIV_,Beta_4,double}
   *
   * The deformation coefficient for the spherical harmonic Y_4_0. This key can
   * be omitted for automatic deformation (`Automatic: true`).
   */
  /**
   * \see_key{key_MC_PT_deformed_betaIV_}
   */
  inline static const Key<double> modi_collider_projectile_deformed_beta4{
      {"Modi", "Collider", "Projectile", "Deformed", "Beta_4"}, {"1.0"}};
  /**
   * \see_key{key_MC_PT_deformed_betaIV_}
   */
  inline static const Key<double> modi_collider_target_deformed_beta4{
      {"Modi", "Collider", "Target", "Deformed", "Beta_4"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_projectile_and_target_
   * \par Defining orientation
   *
   * In the `Orientation` section it is possible to specify the orientation of
   * the nucleus by rotations which are performed about the axes of a coordinate
   * system that is fixed with respect to the nucleus and whose axes are
   * parallel to those of the computational frame before the first rotation.
   * Note that the nucleus is first rotated by phi and then by theta.
   *
   * \optional_key_no_line{key_MC_PT_deformed_orientation_phi_,Phi,double,0.0}
   *
   * The angle by which to rotate the nucleus about the z-axis.
   */
  /**
   * \see_key{key_MC_PT_deformed_orientation_phi_}
   */
  inline static const Key<double>
      modi_collider_projectile_deformed_orientation_phi{
          {"Modi", "Collider", "Projectile", "Deformed", "Orientation", "Phi"},
          0.0,
          {"1.0"}};
  /**
   * \see_key{key_MC_PT_deformed_orientation_phi_}
   */
  inline static const Key<double> modi_collider_target_deformed_orientation_phi{
      {"Modi", "Collider", "Target", "Deformed", "Orientation", "Phi"},
      0.0,
      {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_projectile_and_target_
   * \optional_key_no_line{key_MC_PT_deformed_orientation_theta_,Theta,double,π/2}
   *
   * The angle by which to rotate the nucleus about the rotated x-axis.
   */
  /**
   * \see_key{key_MC_PT_deformed_orientation_theta_}
   */
  inline static const Key<double>
      modi_collider_projectile_deformed_orientation_theta{
          {"Modi", "Collider", "Projectile", "Deformed", "Orientation",
           "Theta"},
          M_PI / 2,
          {"1.0"}};
  /**
   * \see_key{key_MC_PT_deformed_orientation_theta_}
   */
  inline static const Key<double>
      modi_collider_target_deformed_orientation_theta{
          {"Modi", "Collider", "Target", "Deformed", "Orientation", "Theta"},
          M_PI / 2,
          {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_projectile_and_target_
   * \optional_key_no_line{key_MC_PT_deformed_orientation_psi_,Psi,double,0.0}
   *
   * The angle by which to rotate the nucleus about the rotated y-axis.
   */
  /**
   * \see_key{key_MC_PT_deformed_orientation_psi_}
   */
  inline static const Key<double>
      modi_collider_projectile_deformed_orientation_psi{
          {"Modi", "Collider", "Projectile", "Deformed", "Orientation", "Psi"},
          0.0,
          {"1.0"}};
  /**
   * \see_key{key_MC_PT_deformed_orientation_psi_}
   */
  inline static const Key<double> modi_collider_target_deformed_orientation_psi{
      {"Modi", "Collider", "Target", "Deformed", "Orientation", "Psi"},
      0.0,
      {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_projectile_and_target_
   * \optional_key_no_line{key_MC_PT_deformed_orientation_random_,Random_Rotation,bool,false}
   *
   * Whether the created nucleus object should be randomly rotated in space.
   */
  /**
   * \see_key{key_MC_PT_deformed_orientation_random_}
   */
  inline static const Key<bool>
      modi_collider_projectile_deformed_orientation_randomRotation{
          {"Modi", "Collider", "Projectile", "Deformed", "Orientation",
           "Random_Rotation"},
          false,
          {"1.0"}};
  /**
   * \see_key{key_MC_PT_deformed_orientation_random_}
   */
  inline static const Key<bool>
      modi_collider_target_deformed_orientation_randomRotation{
          {"Modi", "Collider", "Target", "Deformed", "Orientation",
           "Random_Rotation"},
          false,
          {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_impact_parameter_
   * \optional_key{key_MC_impact_rnd_reaction_plane_,Random_Reaction_Plane,bool,false}
   *
   * Rotate the direction of the separation of the two nuclei due to the impact
   * parameter with a uniform random angle in the x-y plane.
   */
  /**
   * \see_key{key_MC_impact_rnd_reaction_plane_}
   */
  inline static const Key<bool> modi_collider_impact_randomReactionPlane{
      {"Modi", "Collider", "Impact", "Random_Reaction_Plane"}, false, {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_impact_parameter_
   * \optional_key{key_MC_impact_value_,Value,double,0.0}
   *
   * Fixed value for the impact parameter <b>in fm</b>.
   * \attention If this value is set, all further `Impact` keys are ignored.
   */
  /**
   * \see_key{key_MC_impact_value_}
   */
  inline static const Key<double> modi_collider_impact_value{
      {"Modi", "Collider", "Impact", "Value"}, 0.0, {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_impact_parameter_
   * \optional_key{key_MC_impact_sample_,Sample,string,"quadratic"}
   *
   * - `"uniform"` &rarr; use uniform sampling of the impact parameter
   *   (uniform in \f$b\f$: \f$dP(b) = db\f$)
   * - `"quadratic"` &rarr; use areal (aka quadratic) input sampling (the
   *   probability of an input parameter range is proportional to the area
   *   corresponding to that range, uniform in \f$b^2\f$:
   *   \f$dP(b) = b\,db\f$).
   * - `"custom"` &rarr; requires `Values` and `Yields` to interpolate the
   *   impact parameter distribution and use rejection sampling.
   */
  /**
   * \see_key{key_MC_impact_sample_}
   */
  inline static const Key<std::string> modi_collider_impact_sample{
      {"Modi", "Collider", "Impact", "Sample"}, "quadratic", {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_impact_parameter_
   * <hr>
   * \par Custom sampling
   * \required_key_no_line{key_MC_impact_values_,Values,list of doubles}
   *
   * Values of the impact parameter <b>in fm</b>, with corresponding `Yields`.
   * Must be same length as `Yields`. This key can be omitted if `Sample` is not
   * set to `"custom"`.
   */
  /**
   * \see_key{key_MC_impact_values_}
   */
  inline static const Key<std::vector<double>> modi_collider_impact_values{
      {"Modi", "Collider", "Impact", "Values"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_impact_parameter_
   * \required_key_no_line{key_MC_impact_yields_,Yields,list of doubles}
   *
   * Values of the particle yields, corresponding to `Values`, i.e. the value
   * of the custom distribution at this value. Must be same length as `Values`.
   * This key can be omitted if `Sample` is not set to `"custom"`.
   */
  /**
   * \see_key{key_MC_impact_sample_}
   */
  inline static const Key<std::vector<double>> modi_collider_impact_yields{
      {"Modi", "Collider", "Impact", "Yields"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_impact_parameter_
   * \optional_key{key_MC_impact_range_,Range,list of two doubles,[0.0\,0.0]}
   *
   * A list of minimal and maximal impact parameters <b>in fm</b> between which
   * \f$b\f$ should be chosen. The order of these is not important.
   */
  /**
   * \see_key{key_MC_impact_range_}
   */
  inline static const Key<std::array<double, 2>> modi_collider_impact_range{
      {"Modi", "Collider", "Impact", "Range"}, {{0.0, 0.0}}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_collider_impact_parameter_
   * \optional_key{key_MC_impact_max_,Max,double,0.0}
   *
   * Like `Range: [0.0, Max]`. Note that if both `Range` and `Max` are
   * specified, `Max` takes precedence.
   */
  /**
   * \see_key{key_MC_impact_max_}
   */
  inline static const Key<double> modi_collider_impactMax{
      {"Modi", "Collider", "Impact", "Max"}, 0.0, {"1.0"}};

  /*!\Userguide
   * \page input_modi_sphere_
   * \required_key_no_line{key_MS_radius_,Radius,double}
   *
   * Radius of the sphere <b>in fm</b>.
   */
  /**
   * \see_key{key_MS_radius_}
   */
  inline static const Key<double> modi_sphere_radius{
      {"Modi", "Sphere", "Radius"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_sphere_
   * \required_key{key_MS_temperature_,Temperature,double}
   *
   * Temperature <b>in GeV</b> to sample momenta in the sphere.
   */
  /**
   * \see_key{key_MS_radius_}
   */
  inline static const Key<double> modi_sphere_temperature{
      {"Modi", "Sphere", "Temperature"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_sphere_
   * \required_key{key_MS_start_time_,Start_Time,double}
   *
   * Starting time of sphere calculation.
   */
  /**
   * \see_key{key_MS_start_time_}
   */
  inline static const Key<double> modi_sphere_startTime{
      {"Modi", "Sphere", "Start_Time"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_sphere_
   * \required_key{key_MS_init_mult_,Init_Multiplicities,map<int\,int>}
   *
   * Initial multiplicities per particle species. The value of this key shall be
   * a map of PDG number and amount corresponding to it. Use this key to specify
   * how many particles of each species will be initialized. This key can be
   * omitted if <tt>\ref key_MS_use_thermal_mult_
   * "Use_Thermal_Multiplicities"</tt> is `true`.
   */
  /**
   * \see_key{key_MS_init_mult_}
   */
  inline static const Key<std::map<PdgCode, int>>
      modi_sphere_initialMultiplicities{
          {"Modi", "Sphere", "Init_Multiplicities"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_sphere_
   * \optional_key{key_MS_use_thermal_mult_,Use_Thermal_Multiplicities,bool,false}
   *
   * If this option is set to `true` then <tt>\ref key_MS_init_mult_
   * "Init_Multiplicities"</tt> are ignored and the system is initialized with
   * all particle species of the particle table that belong to the hadron gas
   * equation of state, see HadronGasEos::is_eos_particle(). The multiplicities
   * are sampled from Poisson distributions \f$\mathrm{Poi}(n_i V)\f$, where
   * \f$n_i\f$ are the grand-canonical thermal densities of the corresponding
   * species and \f$V\f$ is the system volume. This option simulates the
   * grand-canonical ensemble, where the number of particles is not fixed from
   * event to event.
   */
  /**
   * \see_key{key_MS_use_thermal_mult_}
   */
  inline static const Key<bool> modi_sphere_useThermalMultiplicities{
      {"Modi", "Sphere", "Use_Thermal_Multiplicities"}, false, {"1.0"}};

  /*!\Userguide
   * \page input_modi_sphere_
   * \optional_key{key_MS_use_bar_chem_pot_,Baryon_Chemical_Potential,double,0.0}
   *
   * Baryon chemical potential \f$\mu_B\f$. This key is used to compute thermal
   * densities \f$n_i\f$ only if <tt>\ref key_MS_use_thermal_mult_
   * "Use_Thermal_Multiplicities"</tt> is `true`.
   */
  /**
   * \see_key{key_MS_use_bar_chem_pot_}
   */
  inline static const Key<double> modi_sphere_baryonChemicalPotential{
      {"Modi", "Sphere", "Baryon_Chemical_Potential"}, 0.0, {"1.0"}};

  /*!\Userguide
   * \page input_modi_sphere_
   * \optional_key{key_MS_strange_chem_pot_,Strange_Chemical_Potential,double,0.0}
   *
   * Strangeness chemical potential \f$\mu_S\f$. This key is used to compute
   * thermal densities \f$n_i\f$ only if <tt>\ref key_MS_use_thermal_mult_
   * "Use_Thermal_Multiplicities"</tt> is `true`.
   */
  /**
   * \see_key{key_MS_strange_chem_pot_}
   */
  inline static const Key<double> modi_sphere_strangeChemicalPotential{
      {"Modi", "Sphere", "Strange_Chemical_Potential"}, 0.0, {"1.0"}};

  /*!\Userguide
   * \page input_modi_sphere_
   * \optional_key{key_MS_charge_chem_pot_,Charge_Chemical_Potential,double,0.0}
   *
   * Charge chemical potential \f$\mu_Q\f$. This key is used to compute
   * thermal densities \f$n_i\f$ only if <tt>\ref key_MS_use_thermal_mult_
   * "Use_Thermal_Multiplicities"</tt> is `true`.
   */
  /**
   * \see_key{key_MS_charge_chem_pot_}
   */
  inline static const Key<double> modi_sphere_chargeChemicalPotential{
      {"Modi", "Sphere", "Charge_Chemical_Potential"}, 0.0, {"1.0"}};

  /*!\Userguide
   * \page input_modi_sphere_
   * \optional_key{key_MS_account_res_widths_,Account_Resonance_Widths,bool,true}
   *
   * This key is considered only in case of thermal initialization and the
   * following two behaviors can be choosen:
   * - `true` &rarr; Account for resonance spectral functions, while computing
   *   multiplicities and sampling masses.
   * - `false` &rarr; Simply use pole masses.
   */
  /**
   * \see_key{key_MS_account_res_widths_}
   */
  inline static const Key<bool> modi_sphere_accountResonanceWidths{
      {"Modi", "Sphere", "Account_Resonance_Widths"}, true, {"1.0"}};

  /*!\Userguide
   * \page input_modi_sphere_
   * \optional_key{key_MS_initial_cond_,Initial_Condition,string,"thermal
   * momenta"}
   *
   * Initial distribution to use for momenta of particles. Mainly used in the
   * expanding universe scenario, options are:
   * - `"thermal momenta"` &rarr; equilibrium Boltzmann distribution
   * - `"IC_ES"` &rarr; off-equilibrium distribution
   * - `"IC_1M"` &rarr; off-equilibrium distribution
   * - `"IC_2M"` &rarr; off-equilibrium distribution
   * - `"IC_Massive"` &rarr; off-equilibrium distribution
   *
   * See \iref{Bazow:2016oky} and \iref{Tindall:2016try} for further
   * explanations about the different distribution functions.
   */
  /**
   * \see_key{key_MS_initial_cond_}
   */
  inline static const Key<std::string> modi_sphere_initialCondition{
      {"Modi", "Sphere", "Initial_Condition"}, "thermal momenta", {"1.0"}};

  /*!\Userguide
   * \page input_modi_sphere_
   * <hr>
   * #### Specifying jets
   *
   * The `Jet` section within the `Sphere` one is used to put a single high
   * energy particle (a "jet") in the center of the system, on an outbound
   * trajectory along the x-axis; if no PDG is specified no jet is produced.
   *
   * \optional_key_no_line{key_MS_jet_jet_pdg_,Jet_PDG,int,no jet}
   *
   * The type of particle to be used as a jet, as given by its PDG code;
   * if none is provided no jet is initialized.
   */
  /**
   * \see_key{key_MS_jet_jet_pdg_}
   */
  inline static const Key<PdgCode> modi_sphere_jet_jetPdg{
      {"Modi", "Sphere", "Jet", "Jet_PDG"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_sphere_
   * \optional_key_no_line{key_MS_jet_jet_momentum_,Jet_Momentum,double,20.0}
   *
   * The initial momentum <b>in GeV</b> to give to the jet particle.
   */
  /**
   * \see_key{key_MS_jet_jet_momentum_}
   */
  inline static const Key<double> modi_sphere_jet_jetMomentum{
      {"Modi", "Sphere", "Jet", "Jet_Momentum"}, 20.0, {"1.0"}};

  /*!\Userguide
   * \page input_modi_box_
   * \required_key_no_line{key_MB_initial_condition_,Initial_Condition,string}
   *
   * Controls initial momentum distribution of particles.
   * - `"peaked momenta"` &rarr; All particles have momentum \f$p=3\,T\f$,
   *   where \f$T\f$ is the temperature. Directions of momenta are uniformly
   *   distributed.
   * - `"thermal momenta"` &rarr; Momenta are sampled from a Maxwell-Boltzmann
   *   distribution.
   */
  /**
   * \see_key{key_MB_initial_condition_}
   */
  inline static const Key<std::string> modi_box_initialCondition{
      {"Modi", "Box", "Initial_Condition"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_box_
   * \required_key{key_MB_length_,Length,double}
   *
   * Length of the cube's edge <b>in fm</b>.
   */
  /**
   * \see_key{key_MB_length_}
   */
  inline static const Key<double> modi_box_length{{"Modi", "Box", "Length"},
                                                  {"1.0"}};

  /*!\Userguide
   * \page input_modi_box_
   * \required_key{key_MB_temperature_,Temperature,double}
   *
   * Temperature <b>in GeV</b> of the box.
   */
  /**
   * \see_key{key_MB_temperature_}
   */
  inline static const Key<double> modi_box_temperature{
      {"Modi", "Box", "Temperature"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_box_
   * \required_key{key_MB_start_time_,Start_Time,double}
   *
   * Starting time of the simulation. All particles in the box are initialized
   * with \f$x^0=\f$`Start_Time`.
   */
  /**
   * \see_key{key_MB_start_time_}
   */
  inline static const Key<double> modi_box_startTime{
      {"Modi", "Box", "Start_Time"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_box_
   * \optional_key{key_MB_equilibration_time_,Equilibration_Time,double, -1.0}
   *
   * Time after which the output of the box is written out. The first time
   * however will be printed. This is useful if one wants to simulate boxes for
   * very long times and knows at which time the box reaches its thermal and
   * chemical equilibrium.
   * The default set to -1 is meaning that output is written from beginning on,
   * if this key is not given.
   */
  /**
   * \see_key{key_MB_equilibration_time_}
   */
  inline static const Key<double> modi_box_equilibrationTime{
      {"Modi", "Box", "Equilibration_Time"}, -1.0, {"1.0"}};

  /*!\Userguide
   * \page input_modi_box_
   * \required_key{key_MB_init_mult_,Init_Multiplicities,map<int\,int>}
   *
   * See &nbsp;
   * <tt>\ref key_MS_init_mult_ "Sphere: Init_Multiplicities"</tt>.
   */
  /**
   * \see_key{key_MB_init_mult_}
   */
  inline static const Key<std::map<PdgCode, int>>
      modi_box_initialMultiplicities{{"Modi", "Box", "Init_Multiplicities"},
                                     {"1.0"}};

  /*!\Userguide
   * \page input_modi_box_
   * \optional_key{key_MB_use_thermal_mult_,Use_Thermal_Multiplicities,bool,false}
   *
   * See &nbsp;
   * <tt>\ref key_MS_use_thermal_mult_
   * "Sphere: Use_Thermal_Multiplicities"</tt>.
   */
  /**
   * \see_key{key_MB_use_thermal_mult_}
   */
  inline static const Key<bool> modi_box_useThermalMultiplicities{
      {"Modi", "Box", "Use_Thermal_Multiplicities"}, false, {"1.0"}};

  /*!\Userguide
   * \page input_modi_box_
   * \optional_key{key_MB_use_bar_chem_pot_,Baryon_Chemical_Potential,double,0.0}
   *
   * See &nbsp;
   * <tt>\ref key_MS_use_bar_chem_pot_ "Sphere: Baryon_Chemical_Potential"</tt>.
   */
  /**
   * \see_key{key_MB_use_bar_chem_pot_}
   */
  inline static const Key<double> modi_box_baryonChemicalPotential{
      {"Modi", "Box", "Baryon_Chemical_Potential"}, 0.0, {"1.0"}};

  /*!\Userguide
   * \page input_modi_box_
   * \optional_key{key_MB_strange_chem_pot_,Strange_Chemical_Potential,double,0.0}
   *
   * See &nbsp;
   * <tt>\ref key_MS_strange_chem_pot_
   * "Sphere: Strange_Chemical_Potential"</tt>.
   */
  /**
   * \see_key{key_MB_strange_chem_pot_}
   */
  inline static const Key<double> modi_box_strangeChemicalPotential{
      {"Modi", "Box", "Strange_Chemical_Potential"}, 0.0, {"1.0"}};

  /*!\Userguide
   * \page input_modi_box_
   * \optional_key{key_MB_charge_chem_pot_,Charge_Chemical_Potential,bool,false}
   *
   * See &nbsp;
   * <tt>\ref key_MS_charge_chem_pot_ "Sphere: Charge_Chemical_Potential"</tt>.
   */
  /**
   * \see_key{key_MB_charge_chem_pot_}
   */
  inline static const Key<double> modi_box_chargeChemicalPotential{
      {"Modi", "Box", "Charge_Chemical_Potential"}, 0.0, {"1.0"}};

  /*!\Userguide
   * \page input_modi_box_
   * \optional_key{key_MB_account_res_widths_,Account_Resonance_Widths,bool,true}
   *
   * See &nbsp;
   * <tt>\ref key_MS_account_res_widths_
   * "Sphere: Account_Resonance_Widths"</tt>.
   *
   * \note
   * Normally, one wants this option `true`. For example, for the detailed
   * balance studies, it is better to account for spectral functions, because
   * then at \f$t=0\f$ one has exactly the expected thermal grand-canonical
   * multiplicities, that can be compared to final ones.  However, by toggling
   * `true` to `false` one can observe the effect of spectral function on the
   * multiplicity. This is useful for understanding the implications of
   * different ways of sampling resonances in hydrodynamics.
   */
  /**
   * \see_key{key_MB_account_res_widths_}
   */
  inline static const Key<bool> modi_box_accountResonanceWidths{
      {"Modi", "Box", "Account_Resonance_Widths"}, true, {"1.0"}};

  /*!\Userguide
   * \page input_modi_box_
   * <hr>
   * #### Specifying jets
   *
   * The `Jet` section can be specified in the `Box` section with the same
   * meaning it has for the `Sphere` modus. It is namely possible to put a
   * jet in the center of the box, on a outbound trajectory along the x-axis.
   *
   * \optional_key_no_line{key_MB_jet_jet_pdg_,Jet_PDG,int,no jet}
   *
   * See &nbsp;
   * <tt>\ref key_MS_jet_jet_pdg_ "Sphere: Jet: Jet_PDG"</tt>.
   */
  /**
   * \see_key{key_MB_jet_jet_pdg_}
   */
  inline static const Key<PdgCode> modi_box_jet_jetPdg{
      {"Modi", "Box", "Jet", "Jet_PDG"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_box_
   * \optional_key_no_line{key_MB_jet_jet_momentum_,Jet_Momentum,double,20.0}
   *
   * See &nbsp;
   * <tt>\ref key_MS_jet_jet_momentum_ "Sphere: Jet: Jet_Momentum"</tt>.
   */
  /**
   * \see_key{key_MB_jet_jet_momentum_}
   */
  inline static const Key<double> modi_box_jet_jetMomentum{
      {"Modi", "Box", "Jet", "Jet_Momentum"}, 20.0, {"1.0"}};

  /*!\Userguide
   * \page input_modi_list_
   * \required_key{key_ML_file_dir_,File_Directory,string}
   *
   * Directory for the external particle lists. Although relative paths to the
   * execution directory should work, you are encouraged to <b>prefer absolute
   * paths</b>.
   */
  /**
   * \see_key{key_ML_file_dir_}
   */
  inline static const Key<std::string> modi_list_fileDirectory{
      {"Modi", "List", "File_Directory"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_list_
   * \required_key{key_ML_file_prefix_,File_Prefix,string}
   *
   * Prefix for the external particle lists file.
   */
  /**
   * \see_key{key_ML_file_prefix_}
   */
  inline static const Key<std::string> modi_list_filePrefix{
      {"Modi", "List", "File_Prefix"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_list_
   * \required_key{key_ML_shift_id_,Shift_Id,int}
   *
   * Starting index for the particle list file(s). To be used to indicate which
   * is the first file to is read.
   */
  /**
   * \see_key{key_ML_shift_id_}
   */
  inline static const Key<std::string> modi_list_shiftId{
      {"Modi", "List", "Shift_Id"}, {"1.0"}};

  /*!\Userguide
   * \page input_modi_listbox_
   * \required_key{key_MLB_length_,Length,double}
   *
   * See &nbsp;
   * <tt>\ref key_MB_length_ "Box: Length"</tt>.
   */
  /**
   * \see_key{key_MB_length_}
   */
  inline static const Key<double> modi_listBox_length{
      {"Modi", "ListBox", "Length"}, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   *
   * ## General output configuration parameters
   *
   * \optional_key_no_line{key_output_out_interval_,Output_Interval,double,
   * \ref key_gen_end_time_ "End_Time"}
   *
   * Defines the period of intermediate output of the status of the simulated
   * system in Standard Output and other output formats which support this
   * functionality.
   */
  /**
   * \see_key{key_output_out_interval_}
   */
  inline static const Key<double> output_outputInterval{
      {"Output", "Output_Interval"}, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * \optional_key{key_output_out_times_,Output_Times,list of doubles,
   * use \ref key_output_out_interval_ "Output_Interval"}
   *
   * Explicitly defines the times where output is generated in the form of
   * a list. This cannot be used in combination with `Output_Interval`. Output
   * times outside the simulation time are ignored and both the initial and
   * final time are always considered. The following example will produce output
   * at event start, event end and at the specified times as long as they are
   * within the simulation time.
   *\verbatim
   Output:
       Output_Times: [-0.1, 0.0, 1.0, 2.0, 10.0]
   \endverbatim
   */
  /**
   * \see_key{key_output_out_times_}
   */
  inline static const Key<std::vector<double>> output_outputTimes{
      {"Output", "Output_Times"}, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * \optional_key{key_output_density_type_,Density_Type,string,"none"}
   *
   * Determines which kind of density is printed into the headers of the
   * collision files. Possible valuesare:
   * - `"hadron"` &rarr; Total hadronic density
   * - `"baryon"` &rarr; Net baryon density
   * - `"baryonic isospin"` &rarr; Baryonic isospin density
   * - `"pion"` &rarr; Pion density
   * - `"none"` &rarr; Do not calculate density, print 0.0
   */
  /**
   * \see_key{key_output_density_type_}
   */
  inline static const Key<std::string> output_densityType{
      {"Output", "Density_Type"}, "none", {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * <hr>
   * ## Output format independently of the specific output content
   *
   * A dedicated subsection in the `Output` section exists for every single
   * output content and dedicated options are described further below. Refer to
   * \ref output_contents_ "output contents" for the list of possible
   * contents. Independently of the content, i.e. in every subsection, it is
   * always necessary (i.e. it is probably desired) to provide the format in
   * which the output should be generated.
   *
   * \optional_key_no_line{key_output_content_format_,Format,list of strings,[]}
   *
   * List of formats for writing particular content. Available formats for every
   * content are listed and described \ref output_contents_ "here", while
   * \ref list_of_output_formats "here" all possible output formats are
   * given.
   *
   * \warning If a `Format` list in a content `section` is not given or it is
   * left empty, i.e. `Format = []`, no output for that given content is
   * produced. Furthermore, if a not existing format is given in the formats
   * list, SMASH is giving a non-fatal error and simply ignoring that format.
   */
  /**
   * \see_key{key_output_content_format_}
   */
  inline static const Key<std::vector<std::string>> output_particles_format{
      {"Output", "Particles", "Format"}, {}, {"1.0"}};
  /**
   * \see_key{key_output_content_format_}
   */
  inline static const Key<std::vector<std::string>> output_collisions_format{
      {"Output", "Collisions", "Format"}, {}, {"1.0"}};
  /**
   * \see_key{key_output_content_format_}
   */
  inline static const Key<std::vector<std::string>> output_dileptons_format{
      {"Output", "Dileptons", "Format"}, {}, {"1.0"}};
  /**
   * \see_key{key_output_content_format_}
   */
  inline static const Key<std::vector<std::string>> output_photons_format{
      {"Output", "Photons", "Format"}, {}, {"1.0"}};
  /**
   * \see_key{key_output_content_format_}
   */
  inline static const Key<std::vector<std::string>>
      output_initialConditions_format{
          {"Output", "Initial_Conditions", "Format"}, {}, {"1.0"}};
  /**
   * \see_key{key_output_content_format_}
   */
  inline static const Key<std::vector<std::string>> output_rivet_format{
      {"Output", "Rivet", "Format"}, {}, {"1.0"}};
  /**
   * \see_key{key_output_content_format_}
   */
  inline static const Key<std::vector<std::string>> output_coulomb_format{
      {"Output", "Coulomb", "Format"}, {}, {"1.0"}};
  /**
   * \see_key{key_output_content_format_}
   */
  inline static const Key<std::vector<std::string>>
      output_thermodynamics_format{
          {"Output", "Thermodynamics", "Format"}, {}, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * <hr>
   * ## Content-specific output options
   * \anchor output_content_specific_options_
   *
   * Every possible content-specific section is documented in the following.
   * Refer to \ref configuring_output_ "this page" for concrete output
   * configuration examples.
   *
   * <hr>
   * ### &diams; Particles
   *
   * \optional_key_no_line{key_output_particles_extended_,Extended,bool,false}
   *
   * &rArr; Incompatible with `Oscar1999`, `VTK`, `HepMC_asciiv3` and
   * `HepMC_treeroot` formats.
   * - `true` &rarr; Print extended information for each particle
   * - `false` &rarr; Regular output for each particle
   */
  /**
   * \see_key{key_output_particles_extended_}
   */
  inline static const Key<bool> output_particles_extended{
      {"Output", "Particles", "Extended"}, false, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * \optional_key_no_line{key_output_particles_only_final_,Only_Final,string,"Yes"}
   *
   * &rArr; Incompatible with `VTK`, `HepMC_asciiv3` and `HepMC_treeroot`
   * - `"Yes"` &rarr; Print only final particle list.
   * - `"IfNotEmpty"` &rarr; Print only final particle list, but only if event
   *   is not empty (i.e. any collisions happened between projectile and
   *   target). Useful to save disk space.
   * - `"No"` &rarr; Particle list at output interval including initial time.
   */
  /**
   * \see_key{key_output_particles_only_final_}
   */
  inline static const Key<std::string> output_particles_onlyFinal{
      {"Output", "Particles", "Only_Final"}, "Yes", {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * <hr>
   * ### &diams; Collisions
   * &rArr; Format `VTK` not available
   *
   * \optional_key_no_line{key_output_collisions_extended_,Extended,bool,false}
   *
   * &rArr; Incompatible with `Oscar1999`, `HepMC_asciiv3` and `HepMC_treeroot`
   * formats.
   * - `true` &rarr; Print extended information for each particle
   * - `false` &rarr; Regular output for each particle
   */
  /**
   * \see_key{key_output_collisions_extended_}
   */
  inline static const Key<bool> output_collisions_extended{
      {"Output", "Collisions", "Extended"}, false, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * \optional_key_no_line{key_output_collisions_print_start_end_,Print_Start_End,bool,false}
   *
   * &rArr; Incompatible with `Root`, `HepMC_asciiv3` and `HepMC_treeroot`
   * formats.
   * - `true` &rarr; Initial and final particle list is printed out
   * - `false` &rarr; Initial and final particle list is not printed out
   */
  /**
   * \see_key{key_output_collisions_print_start_end_}
   */
  inline static const Key<bool> output_collisions_printStartEnd{
      {"Output", "Collisions", "Print_Start_End"}, false, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * <hr>
   * ### &diams; Dileptons
   * &rArr; Only `Oscar1999`, `Oscar2013` and `Binary` formats.
   *
   * \optional_key_no_line{key_output_dileptons_extended_,Extended,bool,false}
   *
   * &rArr; Incompatible with `Oscar1999` format.
   * - `true` &rarr; Print extended information for each particle
   * - `false` &rarr; Regular output for each particle
   */
  /**
   * \see_key{key_output_dileptons_extended_}
   */
  inline static const Key<bool> output_dileptons_extended{
      {"Output", "Dileptons", "Extended"}, false, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * <hr>
   * ### &diams; Photons
   * &rArr; Only `Oscar1999`, `Oscar2013` and `Binary` formats.
   *
   * \optional_key_no_line{key_output_photons_extended_,Extended,bool,false}
   *
   * &rArr; Incompatible with `Oscar1999` format.
   * - `true` &rarr; Print extended information for each particle
   * - `false` &rarr; Regular output for each particle
   */
  /**
   * \see_key{key_output_photons_extended_}
   */
  inline static const Key<bool> output_photons_extended{
      {"Output", "Photons", "Extended"}, false, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * <hr>
   * ### &diams; Initial_Conditions
   * &rArr; Only `Oscar1999`, `Oscar2013`, `Binary`, `ROOT` and `ASCII` (special
   * ASCII IC, see \ref IC_output_user_guide_) formats.
   *
   * \optional_key_no_line{key_output_IC_proper_time_,Proper_Time,double,
   * </tt>\f$f(t_{np})\f$<tt>}
   *
   * Proper time at which hypersurface is created. Its default value depends on
   * the nuclei passing time \f$t_{np}\f$ as follows,
   * \f[
   * f(t_{np})=\begin{cases}
   * \mathrm{\texttt{Lower_Bound}}  & t_{np} \le \mathrm{\texttt{Lower_Bound}}\\
   * t_{np} & t_{np} > \mathrm{\texttt{Lower_Bound}}
   * \end{cases}\;.
   * \f]
   */
  /**
   * \see_key{key_output_IC_proper_time_}
   */
  inline static const Key<double> output_initialConditions_properTime{
      {"Output", "Initial_Conditions", "Proper_Time"}, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * \optional_key_no_line{key_output_IC_lower_bound_,Lower_Bound,double,0.5}
   *
   * Lower bound <b>in fm</b> for the IC proper time if
   * <tt>\ref key_output_IC_proper_time_ "Proper_Time"</tt> is not provided.
   */
  /**
   * \see_key{key_output_IC_lower_bound_}
   */
  inline static const Key<double> output_initialConditions_lowerBound{
      {"Output", "Initial_Conditions", "Lower_Bound"}, 0.5, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * \optional_key_no_line{key_output_IC_rapidity_cut_,Rapidity_Cut,double,
   * </tt>No cut is done<tt>}
   *
   * If set, employ a rapidity cut for particles contributing to the initial
   * conditions for hydrodynamics. A positive value is expected and the cut is
   * employed symmetrically around 0. Only particles characterized by
   * \f$|\mathrm{\texttt{Rapidity_Cut}}|<y\f$ are printed to the
   * output file.
   */
  /**
   * \see_key{key_output_IC_rapidity_cut_}
   */
  inline static const Key<double> output_initialConditions_rapidityCut{
      {"Output", "Initial_Conditions", "Rapidity_Cut"}, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * \optional_key_no_line{key_output_IC_pt_cut_,pT_Cut,double,
   * </tt>No cut is done<tt>}
   *
   * If set, employ a transverse momentum cut for particles contributing to the
   * initial conditions for hydrodynamics. A positive value is expected. Only
   * particles characterized by \f$0<p_T<\mathrm{\texttt{pT_Cut}}\f$ are printed
   * to the output file.
   */
  /**
   * \see_key{key_output_IC_pt_cut_}
   */
  inline static const Key<double> output_initialConditions_pTCut{
      {"Output", "Initial_Conditions", "pT_Cut"}, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * \optional_key_no_line{key_output_IC_extended_,Extended,bool,false}
   *
   * &rArr; Incompatible with `Oscar1999`, `ROOT` and `ASCII` formats.
   * - `true` &rarr' Print extended information for each particle
   * - `false` &rarr' Regular output for each particle
   */
  /**
   * \see_key{key_output_IC_extended_}
   */
  inline static const Key<bool> output_initialConditions_extended{
      {"Output", "Initial_Conditions", "Extended"}, false, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * <hr> \anchor input_output_rivet_
   * ### &diams; Rivet
   * &rArr; Only `YODA` format (see \ref rivet_output_user_guide_ "here" for
   * more information about the format).
   *
   * \note In the following, <b>no default</b> means that, if the key is
   *       omitted, Rivet default behavior will be used.
   *
   * \optional_key_no_line{key_output_rivet_paths_,Paths,list of strings,
   * </tt><b>no default</b><tt>}
   *
   * This key specifies the directories that Rivet will search for analyses
   * and data files related to the analyses.
   */
  /**
   * \see_key{key_output_rivet_paths_}
   */
  inline static const Key<std::vector<std::string>> output_rivet_paths{
      {"Output", "Rivet", "Paths"}, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * \optional_key_no_line{key_output_rivet_analyses_,Analyses,list of strings,
   * </tt><b>no default</b><tt>}
   *
   * This key specifies the analyses (including possible options) to add to the
   * Rivet analysis.
   */
  /**
   * \see_key{key_output_rivet_analyses_}
   */
  inline static const Key<std::vector<std::string>> output_rivet_analyses{
      {"Output", "Rivet", "Analyses"}, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * \optional_key_no_line{key_output_rivet_preloads_,Preloads,list of strings,
   * </tt><b>no default</b><tt>}
   *
   * Specify data files to read into Rivet (e.g., centrality calibrations) at
   * start-up.
   */
  /**
   * \see_key{key_output_rivet_preloads_}
   */
  inline static const Key<std::vector<std::string>> output_rivet_preloads{
      {"Output", "Rivet", "Preloads"}, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * \optional_key_no_line{key_output_rivet_logging_,Logging,map<string\,string>,
   * </tt><b>no default</b><tt>}
   *
   * Specifies log levels for various parts of Rivet, including analyses. Each
   * entry is a log name followed by a log level (one among `"TRACE"`,
   * `"DEBUG"`, `"INFO"`, `"WARN"`, `"ERROR"`, and `"FATAL"`).
   */
  /**
   * \see_key{key_output_rivet_logging_}
   */
  inline static const Key<std::map<std::string, std::string>>
      output_rivet_logging{{"Output", "Rivet", "Logging"}, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * \optional_key_no_line{key_output_rivet_ignore_beams_,Ignore_Beams,bool,true}
   *
   * Ask Rivet to not validate beams before running analyses. This is needed if
   * you use the <tt>\ref key_MC_fermi_motion_ "Fermi_Motion"</tt> option that
   * disrupts the collision energy event-by-event.
   */
  /**
   * \see_key{key_output_rivet_ignore_beams_}
   */
  inline static const Key<bool> output_rivet_ignoreBeams{
      {"Output", "Rivet", "Ignore_Beams"}, true, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * \optional_key_no_line{key_output_rivet_cross_sections_,Cross_Section,
   * list of two doubles,</tt><b>no default</b><tt>}
   *
   * Set the cross-section <b>in pico-barns</b>.
   */
  /**
   * \see_key{key_output_rivet_cross_sections_}
   */
  inline static const Key<std::array<double, 2>> output_rivet_crossSection{
      {"Output", "Rivet", "Cross_Section"}, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   *
   * #### Weights keys
   *
   * Some operations about weights can be customized in the `Weights` section.
   *
   * \optional_key_no_line{key_output_rivet_weights_no_multi_,No_Multi,bool,
   * </tt><b>no default</b><tt>}
   *
   * Ask Rivet not to do multi-weight processing.
   */
  /**
   * \see_key{key_output_rivet_weights_no_multi_}
   */
  inline static const Key<std::array<double, 2>> output_rivet_weights_noMulti{
      {"Output", "Rivet", "Weights", "No_Multi"}, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * \optional_key_no_line{key_output_rivet_weights_nominal_,Nominal,string,
   * </tt><b>no default</b><tt>}
   *
   * The nominal weight name.
   */
  /**
   * \see_key{key_output_rivet_weights_nominal_}
   */
  inline static const Key<std::string> output_rivet_weights_nominal{
      {"Output", "Rivet", "Weights", "Nominal"}, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * \optional_key_no_line{key_output_rivet_weights_select_,Select,
   * list of strings, </tt><b>no default</b><tt>}
   *
   * Select these weights for processing.
   */
  /**
   * \see_key{key_output_rivet_weights_select_}
   */
  inline static const Key<std::vector<std::string>> output_rivet_weights_select{
      {"Output", "Rivet", "Weights", "Select"}, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * \optional_key_no_line{key_output_rivet_weights_deselect_,Deselect,
   * list of strings, </tt><b>no default</b><tt>}
   *
   * De-select these weights for processing.
   */
  /**
   * \see_key{key_output_rivet_weights_deselect_}
   */
  inline static const Key<std::vector<std::string>>
      output_rivet_weights_deselect{{"Output", "Rivet", "Weights", "Deselect"},
                                    {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * \optional_key_no_line{key_output_rivet_weights_nlo_smearing_,NLO_Smearing,
   * double, </tt><b>no default</b><tt>}
   *
   * Smearing histogram binning by given fraction of bin widths to avoid NLO
   * counter events to flow into neighboring bin.
   */
  /**
   * \see_key{key_output_rivet_weights_nlo_smearing_}
   */
  inline static const Key<double> output_rivet_weights_nloSmearing{
      {"Output", "Rivet", "Weights", "NLO_Smearing"}, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * \optional_key_no_line{key_output_rivet_weights_cap_,Cap,double,
   * </tt><b>no default</b><tt>}
   *
   * Cap weights to this value.
   */
  /**
   * \see_key{key_output_rivet_weights_cap_}
   */
  inline static const Key<double> output_rivet_weights_cap{
      {"Output", "Rivet", "Weights", "Cap"}, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * <hr>
   * ### &diams; Coulomb
   * &rArr; Only `VTK` format.
   *
   * No content-specific output options, apart from the <tt>\ref
   * key_output_content_format_ "Format"</tt> key which accept `["VTK"]` value
   * only.
   */

  /*!\Userguide
   * \page input_output_
   * <hr> \anchor input_output_thermodynamics_
   * ### &diams; Thermodynamics
   *
   * The user can print thermodynamical quantities
   * -# on the spatial lattice to VTK output;
   * -# on the spatial lattice to ASCII output;
   * -# at a given point to ASCII output;
   * -# averaged over all particles to ASCII output.
   *
   * <b>About 1 and 2:</b> Note that this output requires a lattice, which needs
   * to be enabled in the conguration file and is regulated by the options of
   * \ref input_lattice_. See \ref output_vtk_lattice_ for further information.
   *
   * <b>About 3 and 4:</b> See \ref thermodyn_output_user_guide_ for further
   * information.
   *
   * \optional_key_no_line{key_output_thermo_type_,Type,string,"baryon"}
   *
   * Particle type taken into consideration, one among
   * - `"hadron"`
   * - `"baryon"` (corresponds to "net baryon")
   * - `"baryonic isospin"`
   * - `"pion"`
   * - `"none"`
   * - `"total isospin"`
   */
  /**
   * \see_key{key_output_thermo_type_}
   */
  inline static const Key<DensityType> output_thermodynamics_type{
      {"Output", "Thermodynamics", "Type"}, DensityType::Baryon, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * \optional_key_no_line{key_output_thermo_quantities_,Quantities,
   * list of strings,[]}
   *
   * List of thermodynamic quantities that are printed to the output.
   * Possible quantities are:
   * - `"rho_eckart"` &rarr; Eckart rest frame density.
   * - `"tmn"` &rarr; Energy-momentum tensor \f$T^{\mu\nu}(t,x,y,z)\f$.
   * - `"tmn_landau"` &rarr; Energy-momentum tensor in the Landau rest frame.
   *   This tensor is computed by boosting \f$T^{\mu\nu}(t,x,y,z)\f$ to the
   *   local rest frame, where \f$T^{0i}\f$ = 0.
   * - `"landau_velocity"` &rarr; Velocity of the Landau rest frame. The
   *   velocity is obtained from the energy-momentum tensor
   *   \f$T^{\mu\nu}(t,x,y,z)\f$ by solving the generalized eigenvalue equation
   *   \f$(T^{\mu\nu} - \lambda g^{\mu\nu})u_{\mu}=0\f$.
   * - `"j_QBS"` &rarr; Electric (Q), baryonic (B) and strange (S) currents
   *   \f$j^{\mu}_{QBS}(t,x,y,z) \f$; note that all currents are given in units
   *   of "number of charges"; multiply the electric current by the elementary
   *   charge \f$\sqrt{4 \pi \alpha_{EM}} \f$ for charge units.
   */
  /**
   * \see_key{key_output_thermo_type_}
   */
  inline static const Key<std::vector<std::string>>
      output_thermodynamics_quantites{
          {"Output", "Thermodynamics", "Quantities"}, {}, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * \optional_key_no_line{key_output_thermo_position_,Position,
   * list of 3 doubles,[0.0\, 0.0\, 0.0]}
   *
   * Point at which thermodynamic quantities are computed.
   */
  /**
   * \see_key{key_output_thermo_position_}
   */
  inline static const Key<std::array<double, 3>> output_thermodynamics_position{
      {"Output", "Thermodynamics", "Position"}, {{0.0, 0.0, 0.0}}, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * \optional_key_no_line{key_output_thermo_smearing_,Smearing,bool,true}
   *
   * Using Gaussian smearing for computing thermodynamic quantities or not. This
   * triggers whether thermodynamic quantities are evaluated at a fixed point
   * (`true`) or summed over all particles (`false`).
   * - `true` &rarr; smearing applied
   * - `false` &rarr; smearing not applied
   *
   * The contribution to the energy-momentum tensor and current (be it electric,
   * baryonic or strange) from a single particle in its rest frame is:
   * \f[\begin{eqnarray}
   * j^{\mu} = B \frac{p_0^{\mu}}{p_0^0} W \\
   * T^{\mu \nu} = \frac{p_0^{\mu}p_0^{\nu}}{p_0^0} W
   * \end{eqnarray}
   * \f]
   * with B being the charge of interest and W being the weight given to this
   * particle. Normally, if one computes thermodynamic quantities at a point,
   * smearing should be applied, and then \f$W\f$ takes on the following shape:
   * \f[
   * W = (2 \pi \sigma^2)^{-3/2} \exp\left(
   * - \frac{(\mathbf{r}-\mathbf{r_0(t)})^2}{2\sigma^2}
   * \right)\f]
   * It can however be useful to compute the thermodynamic quantities of all
   * particles in a box with \f$W=1\f$, which would correspond to <tt>"Smearing:
   * false"</tt>. Note that using this option changes the units of the
   * thermodynamic quantities, as they are no longer spatially normalized. One
   * should divide this quantity by the volume of the box to restore units to
   * the correct ones.
   */
  /**
   * \see_key{key_output_thermo_smearing_}
   */
  inline static const Key<bool> output_thermodynamics_smearing{
      {"Output", "Thermodynamics", "Smearing"}, true, {"1.0"}};

  /*!\Userguide
   * \page input_output_
   * \optional_key_no_line{key_output_thermo_only_part_,Only_Participants,bool,false}
   *
   * If set to `true`, only participants are included in the computation of the
   * energy momentum tensor and of the Eckart currents. In this context, a
   * hadron is considered as a participant if it had at least one collision.
   * When using \ref input_potentials_ "Potentials" this option must be either
   * left unset or set to `false`. The reason behing this limitation is that in
   * this case hadrons can influence the evolution of the system even without
   * collisions.
   */
  /**
   * \see_key{key_output_thermo_only_part_}
   */
  inline static const Key<bool> output_thermodynamics_onlyParticipants{
      {"Output", "Thermodynamics", "Only_Participants"}, false, {"1.0"}};

  /*!\Userguide
   * \page input_lattice_
   * \optional_key{key_lattice_sizes_,Sizes,list of 3 doubles,
   * </tt>depends on <tt>\ref key_gen_modus_ "Modus"}
   *
   * Sizes of lattice in x, y, z directions <b>in fm</b>.
   */
  /**
   * \see_key{key_lattice_sizes_}
   */
  inline static const Key<std::array<double, 3>> lattice_sizes{
      {"Lattice", "Sizes"}, {"1.0"}};

  /*!\Userguide
   * \page input_lattice_
   * \optional_key{key_lattice_cell_number_,Cell_Number,list of 3 ints,
   * </tt>depends on <tt>\ref key_gen_modus_ "Modus"}
   *
   * Number of cells in x, y, z directions.
   */
  /**
   * \see_key{key_lattice_cell_number_}
   */
  inline static const Key<std::array<int, 3>> lattice_cellNumber{
      {"Lattice", "Cell_Number"}, {"1.0"}};

  /*!\Userguide
   * \page input_lattice_
   * \optional_key{key_lattice_origin_,Origin,list of 3 doubles,
   * </tt>depends on <tt>\ref key_gen_modus_ "Modus"}
   *
   * Coordinates of the left, down, near corner of the lattice <b>in fm</b>.
   */
  /**
   * \see_key{key_lattice_origin_}
   */
  inline static const Key<std::array<double, 3>> lattice_origin{
      {"Lattice", "Origin"}, {"1.0"}};

  /*!\Userguide
   * \page input_lattice_
   * \optional_key{key_lattice_periodic_,Periodic,bool,
   * (\ref key_gen_modus_ "Modus" == "Box")}
   *
   * Use periodic continuation or not. With periodic continuation
   * \f$(x,y,z) + (i\cdot L_x,\,j\cdot L_y,\,k\cdot L_z) \equiv (x,y,z)\f$
   * with \f$i,\,j,\,k\in\mathbb{Z}\f$ and \f$L_x,\,L_y,\,L_z\f$ being the
   * lattice sizes.
   */
  /**
   * \see_key{key_lattice_periodic_}
   */
  inline static const Key<bool> lattice_periodic{{"Lattice", "Periodic"},
                                                 {"1.0"}};

  /*!\Userguide
   * \page input_lattice_
   * \optional_key{key_lattice_pot_affect_threshold_,Potentials_Affect_Thresholds,bool,false}
   *
   * Include potential effects, since mean field potentials change the threshold
   * energies of the actions.
   */
  /**
   * \see_key{key_lattice_pot_affect_threshold_}
   */
  inline static const Key<bool> lattice_potentialsAffectThreshold{
      {"Lattice", "Potentials_Affect_Thresholds"}, false, {"1.0"}};

  /*!\Userguide
   * \page input_potentials_skyrme_
   * \required_key{key_potentials_skyrme_a_,Skyrme_A,double}
   *
   * Parameter \f$A\f$ of Skyrme potential <b>in MeV</b>.
   */
  /**
   * \see_key{key_potentials_skyrme_a_}
   */
  inline static const Key<double> potentials_skyrme_skyrmeA{
      {"Potentials", "Skyrme", "Skyrme_A"}, {"1.0"}};

  /*!\Userguide
   * \page input_potentials_skyrme_
   * \required_key{key_potentials_skyrme_b_,Skyrme_B,double}
   *
   * Parameter \f$B\f$ of Skyrme potential <b>in MeV</b>.
   */
  /**
   * \see_key{key_potentials_skyrme_b_}
   */
  inline static const Key<double> potentials_skyrme_skyrmeB{
      {"Potentials", "Skyrme", "Skyrme_B"}, {"1.0"}};

  /*!\Userguide
   * \page input_potentials_skyrme_
   * \required_key{key_potentials_skyrme_tau_,Skyrme_Tau,double}
   *
   * Parameter \f$\tau\f$ of Skyrme potential.
   *
   */
  /**
   * \see_key{key_potentials_skyrme_tau_}
   */
  inline static const Key<double> potentials_skyrme_skyrmeTau{
      {"Potentials", "Skyrme", "Skyrme_Tau"}, {"1.0"}};

  /*!\Userguide
   * \page input_potentials_symmetry_
   * \required_key{key_potentials_symmetry_s_pot_,S_pot,double}
   *
   * Parameter \f$S_{pot}\f$ of symmetry potential <b>in MeV</b>.
   */
  /**
   * \see_key{key_potentials_symmetry_s_pot_}
   */
  inline static const Key<double> potentials_symmetry_sPot{
      {"Potentials", "Symmetry", "S_pot"}, {"1.0"}};

  /*!\Userguide
   * \page input_potentials_symmetry_
   * \optional_key{key_potentials_symmetry_gamma_,gamma,double,
   * </tt>do not consider last term in \f$S(\rho_B)\f$<tt>}
   *
   * Exponent \f$\gamma\f$ in formula for \f$S(\rho_B)\f$. If `gamma` is
   * specified, the baryon density dependence is included in the potential.
   * Otherwise only the first term of the potential will be taken into account.
   */
  /**
   * \see_key{key_potentials_symmetry_gamma_}
   */
  inline static const Key<double> potentials_symmetry_gamma{
      {"Potentials", "Symmetry", "gamma"}, {"1.0"}};

  /*!\Userguide
   * \page input_potentials_VDF_
   * \required_key{key_potentials_vdf_sat_rhoB_,Sat_rhoB,double}
   *
   * The saturation density of nuclear matter <b>in 1/fm³</b>.
   */
  /**
   * \see_key{key_potentials_symmetry_gamma_}
   */
  inline static const Key<double> potentials_vdf_satRhoB{
      {"Potentials", "VDF", "Sat_rhoB"}, {"1.0"}};

  /*!\Userguide
   * \page input_potentials_VDF_
   * \required_key{key_potentials_vdf_coeffs_,Coeffs,list of doubles}
   *
   * Parameters \f$C_i\f$ of the VDF potential <b>in MeV</b>.
   */
  /**
   * \see_key{key_potentials_vdf_coeffs_}
   */
  inline static const Key<std::vector<double>> potentials_vdf_coeffs{
      {"Potentials", "VDF", "Coeffs"}, {"1.0"}};

  /*!\Userguide
   * \page input_potentials_VDF_
   * \required_key{key_potentials_vdf_powers_,Powers,double}
   *
   * Parameters \f$b_i\f$ of the VDF potential.
   *
   * \warning
   * You need to provide as many entries for `Powers` as provided for `Coeffs`.
   */
  /**
   * \see_key{key_potentials_vdf_powers_}
   */
  inline static const Key<std::vector<double>> potentials_vdf_powers{
      {"Potentials", "VDF", "Powers"}, {"1.0"}};

  /*!\Userguide
   * \page input_potentials_coulomb_
   * \required_key{key_potentials_coulomb_r_cut_,R_Cut,double}
   *
   * The value at which the integration volume is cut.
   */
  /**
   * \see_key{key_potentials_coulomb_r_cut_}
   */
  inline static const Key<std::vector<double>> potentials_coulomb_rCut{
      {"Potentials", "Coulomb", "R_Cut"}, {"1.0"}};

  /*!\Userguide
   * \page input_forced_thermalization_
   * \required_key{key_forced_therm_cell_number_,Cell_Number,list of 3 doubles}
   *
   * Number of cells in each direction (x,y,z).
   */
  /**
   * \see_key{key_forced_therm_cell_number_}
   */
  inline static const Key<std::array<int, 3>> forcedThermalization_cellNumber{
      {"Forced_Thermalization", "Cell_Number"}, {"1.0"}};

  /*!\Userguide
   * \page input_forced_thermalization_
   * \required_key{key_forced_therm_critical_edens_,Critical_Edens,double}
   *
   * Critical energy density <b>in GeV/fm³</b> above which forced thermalization
   * is applied.
   */
  /**
   * \see_key{key_forced_therm_critical_edens_}
   */
  inline static const Key<double> forcedThermalization_criticalEDensity{
      {"Forced_Thermalization", "Critical_Edens"}, {"1.0"}};

  /*!\Userguide
   * \page input_forced_thermalization_
   * \required_key{key_forced_therm_start_time_,Start_Time,double}
   *
   * Time <b>in fm/c</b> after which forced thermalization may be applied, if
   * the energy density is sufficiently high.
   */
  /**
   * \see_key{key_forced_therm_start_time_}
   */
  inline static const Key<double> forcedThermalization_startTime{
      {"Forced_Thermalization", "Start_Time"}, {"1.0"}};

  /*!\Userguide
   * \page input_forced_thermalization_
   * \required_key{key_forced_therm_timestep_,Timestep,double}
   *
   * Timestep of thermalization <b>in fm/c</b>.
   */
  /**
   * \see_key{key_forced_therm_timestep_}
   */
  inline static const Key<double> forcedThermalization_timestep{
      {"Forced_Thermalization", "Timestep"}, {"1.0"}};

  /*!\Userguide
   * \page input_forced_thermalization_
   * \optional_key{key_forced_therm_algorithm_,Algorithm,string,"biased BF"}
   *
   * Algorithm applied to enforce thermalization, see
   * \iref{Oliinychenko:2016vkg} for more details.
   * - `"unbiased BF"` &rarr; slowest, but theoretically most robust
   * - `"biased BF"` &rarr; faster, but theoretically less robust
   * - `"mode sampling"` &rarr; fastest, but least robust
   */
  /**
   * \see_key{key_forced_therm_algorithm_}
   */
  inline static const Key<std::string> forcedThermalization_algorithm{
      {"Forced_Thermalization", "Algorithm"}, "biased BF", {"1.0"}};

  /*!\Userguide
   * \page input_forced_thermalization_
   * \optional_key{key_forced_therm_microcanonical_,Microcanonical,bool,false}
   *
   * Enforce energy conservation or not as part of sampling algorithm. Relevant
   * for biased and unbiased Becattini-Ferroni (BF) algorithms. If this option
   * is on, samples with energies deviating too far from the initial one will be
   * rejected. This is different from simple energy and momentum
   * renormalization, which is done in the end anyway. If energy conservation
   * is enforced at sampling, the distributions become microcanonical instead
   * of canonical. One particular effect is that multiplicity distributions
   * become narrower.
   *
   * The downside of having this option on is that the sampling takes
   * significantly longer time.
   */
  /**
   * \see_key{key_forced_therm_microcanonical_}
   */
  inline static const Key<bool> forcedThermalization_microcanonical{
      {"Forced_Thermalization", "Microcanonical"}, false, {"1.0"}};

  /*!\Userguide
   * \page input_forced_thermalization_
   * \required_key{key_forced_therm_lattice_sizes_,Lattice_Sizes,list of 3
   * doubles}
   *
   * The lattice is placed such that the center is [0.0,0.0,0.0].
   * If one wants to have a central cell with center at [0.0,0.0,0.0] then
   * number of cells should be odd (2k+1) in every direction.
   *
   * `Lattice_Sizes` is required for all modi, except the `"Box"` modus. In
   * case of `"Box"` modus, the lattice is set up automatically to match the box
   * size, and the user should not (and is not allowed to) specify it.
   */
  /**
   * \see_key{key_forced_therm_lattice_sizes_}
   */
  inline static const Key<std::array<double, 3>>
      forcedThermalization_latticeSizes{
          {"Forced_Thermalization", "Lattice_Sizes"}, {"1.0"}};

  /// Alias for the type to be used in the list of keys.
  using key_references_variant = std::variant<
      std::reference_wrapper<const Key<bool>>,
      std::reference_wrapper<const Key<int>>,
      std::reference_wrapper<const Key<double>>,
      std::reference_wrapper<const Key<std::string>>,
      std::reference_wrapper<const Key<std::array<int, 3>>>,
      std::reference_wrapper<const Key<std::array<double, 2>>>,
      std::reference_wrapper<const Key<std::array<double, 3>>>,
      std::reference_wrapper<const Key<std::vector<double>>>,
      std::reference_wrapper<const Key<std::vector<std::string>>>,
      std::reference_wrapper<const Key<std::map<PdgCode, int>>>,
      std::reference_wrapper<const Key<std::map<std::string, std::string>>>,
      std::reference_wrapper<const Key<einhard::LogLevel>>,
      std::reference_wrapper<const Key<DensityType>>,
      std::reference_wrapper<const Key<DerivativesMode>>,
      std::reference_wrapper<const Key<ExpansionMode>>,
      std::reference_wrapper<const Key<MultiParticleReactionsBitSet>>,
      std::reference_wrapper<const Key<PdgCode>>,
      std::reference_wrapper<const Key<ReactionsBitSet>>,
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
      std::cref(log_yamlConfiguration),
      std::cref(collTerm_twoToOne),
      std::cref(collTerm_includedTwoToTwo),
      std::cref(collTerm_multiParticleReactions),
      std::cref(collTerm_forceDecaysAtEnd),
      std::cref(collTerm_noCollisions),
      std::cref(collTerm_nnbarTreatment),
      std::cref(collTerm_useAQM),
      std::cref(collTerm_resonanceLifetimeModifier),
      std::cref(collTerm_stringsWithProbability),
      std::cref(collTerm_elasticCrossSection),
      std::cref(collTerm_isotropic),
      std::cref(collTerm_maximumCrossSection),
      std::cref(collTerm_fixedMinCellLength),
      std::cref(collTerm_crossSectionScaling),
      std::cref(collTerm_additionalElasticCrossSection),
      std::cref(collTerm_includeDecaysAtTheEnd),
      std::cref(collTerm_elasticNNCutoffSqrts),
      std::cref(collTerm_strings),
      std::cref(collTerm_collisionCriterion),
      std::cref(collTerm_onlyWarnForHighProbability),
      std::cref(collTerm_pauliBlocking_spatialAveragingRadius),
      std::cref(collTerm_pauliBlocking_gaussianCutoff),
      std::cref(collTerm_pauliBlocking_momentumAveragingRadius),
      std::cref(collTerm_stringParam_stringTension),
      std::cref(collTerm_stringParam_gluonBeta),
      std::cref(collTerm_stringParam_gluonPMin),
      std::cref(collTerm_stringParam_quarkAlpha),
      std::cref(collTerm_stringParam_quarkBeta),
      std::cref(collTerm_stringParam_strangeSuppression),
      std::cref(collTerm_stringParam_diquarkSuppression),
      std::cref(collTerm_stringParam_sigmaPerp),
      std::cref(collTerm_stringParam_stringZA),
      std::cref(collTerm_stringParam_stringZB),
      std::cref(collTerm_stringParam_separateFragmentBaryon),
      std::cref(collTerm_stringParam_stringZALeading),
      std::cref(collTerm_stringParam_stringZBLeading),
      std::cref(collTerm_stringParam_stringSigmaT),
      std::cref(collTerm_stringParam_formTimeFactor),
      std::cref(collTerm_stringParam_powerParticleFormation),
      std::cref(collTerm_stringParam_formationTime),
      std::cref(collTerm_stringParam_mDependentFormationTimes),
      std::cref(collTerm_stringParam_probabilityPToDUU),
      std::cref(collTerm_stringParam_popcornRate),
      std::cref(collTerm_dileptons_decays),
      std::cref(collTerm_photons_fractionalPhotons),
      std::cref(collTerm_photons_twoToTwoScatterings),
      std::cref(collTerm_photons_bremsstrahlung),
      std::cref(modi_collider_sqrtSNN),
      std::cref(modi_collider_eKin),
      std::cref(modi_collider_eTot),
      std::cref(modi_collider_pLab),
      std::cref(modi_collider_calculationFrame),
      std::cref(modi_collider_fermiMotion),
      std::cref(modi_collider_collisionWithinNucleus),
      std::cref(modi_collider_initialDistance),
      std::cref(modi_collider_projectile_particles),
      std::cref(modi_collider_target_particles),
      std::cref(modi_collider_projectile_diffusiveness),
      std::cref(modi_collider_target_diffusiveness),
      std::cref(modi_collider_projectile_radius),
      std::cref(modi_collider_target_radius),
      std::cref(modi_collider_projectile_saturationDensity),
      std::cref(modi_collider_target_saturationDensity),
      std::cref(modi_collider_projectile_eTot),
      std::cref(modi_collider_target_eTot),
      std::cref(modi_collider_projectile_eKin),
      std::cref(modi_collider_target_eKin),
      std::cref(modi_collider_projectile_pLab),
      std::cref(modi_collider_target_pLab),
      std::cref(modi_collider_projectile_custom_fileDirectory),
      std::cref(modi_collider_target_custom_fileDirectory),
      std::cref(modi_collider_projectile_custom_fileName),
      std::cref(modi_collider_target_custom_fileName),
      std::cref(modi_collider_projectile_deformed_automatic),
      std::cref(modi_collider_target_deformed_automatic),
      std::cref(modi_collider_projectile_deformed_beta2),
      std::cref(modi_collider_target_deformed_beta2),
      std::cref(modi_collider_projectile_deformed_beta3),
      std::cref(modi_collider_target_deformed_beta3),
      std::cref(modi_collider_projectile_deformed_beta4),
      std::cref(modi_collider_target_deformed_beta4),
      std::cref(modi_collider_projectile_deformed_orientation_phi),
      std::cref(modi_collider_target_deformed_orientation_phi),
      std::cref(modi_collider_projectile_deformed_orientation_psi),
      std::cref(modi_collider_target_deformed_orientation_psi),
      std::cref(modi_collider_projectile_deformed_orientation_theta),
      std::cref(modi_collider_target_deformed_orientation_theta),
      std::cref(modi_collider_projectile_deformed_orientation_randomRotation),
      std::cref(modi_collider_target_deformed_orientation_randomRotation),
      std::cref(modi_collider_impact_randomReactionPlane),
      std::cref(modi_collider_impact_value),
      std::cref(modi_collider_impact_sample),
      std::cref(modi_collider_impact_values),
      std::cref(modi_collider_impact_yields),
      std::cref(modi_collider_impact_range),
      std::cref(modi_collider_impactMax),
      std::cref(modi_sphere_radius),
      std::cref(modi_sphere_temperature),
      std::cref(modi_sphere_startTime),
      std::cref(modi_sphere_initialMultiplicities),
      std::cref(modi_sphere_useThermalMultiplicities),
      std::cref(modi_sphere_baryonChemicalPotential),
      std::cref(modi_sphere_strangeChemicalPotential),
      std::cref(modi_sphere_chargeChemicalPotential),
      std::cref(modi_sphere_accountResonanceWidths),
      std::cref(modi_sphere_initialCondition),
      std::cref(modi_sphere_jet_jetPdg),
      std::cref(modi_sphere_jet_jetMomentum),
      std::cref(modi_box_initialCondition),
      std::cref(modi_box_length),
      std::cref(modi_box_temperature),
      std::cref(modi_box_startTime),
      std::cref(modi_box_equilibrationTime),
      std::cref(modi_box_initialMultiplicities),
      std::cref(modi_box_useThermalMultiplicities),
      std::cref(modi_box_baryonChemicalPotential),
      std::cref(modi_box_strangeChemicalPotential),
      std::cref(modi_box_chargeChemicalPotential),
      std::cref(modi_box_accountResonanceWidths),
      std::cref(modi_box_jet_jetPdg),
      std::cref(modi_box_jet_jetMomentum),
      std::cref(modi_list_fileDirectory),
      std::cref(modi_list_filePrefix),
      std::cref(modi_list_shiftId),
      std::cref(modi_listBox_length),
      std::cref(output_outputInterval),
      std::cref(output_outputTimes),
      std::cref(output_densityType),
      std::cref(output_particles_format),
      std::cref(output_collisions_format),
      std::cref(output_dileptons_format),
      std::cref(output_photons_format),
      std::cref(output_initialConditions_format),
      std::cref(output_rivet_format),
      std::cref(output_coulomb_format),
      std::cref(output_thermodynamics_format),
      std::cref(output_particles_extended),
      std::cref(output_particles_onlyFinal),
      std::cref(output_collisions_extended),
      std::cref(output_collisions_printStartEnd),
      std::cref(output_dileptons_extended),
      std::cref(output_photons_extended),
      std::cref(output_initialConditions_properTime),
      std::cref(output_initialConditions_lowerBound),
      std::cref(output_initialConditions_rapidityCut),
      std::cref(output_initialConditions_pTCut),
      std::cref(output_initialConditions_extended),
      std::cref(output_rivet_paths),
      std::cref(output_rivet_analyses),
      std::cref(output_rivet_preloads),
      std::cref(output_rivet_preloads),
      std::cref(output_rivet_preloads),
      std::cref(output_rivet_crossSection),
      std::cref(output_rivet_weights_noMulti),
      std::cref(output_rivet_weights_nominal),
      std::cref(output_rivet_weights_select),
      std::cref(output_rivet_weights_deselect),
      std::cref(output_rivet_weights_nloSmearing),
      std::cref(output_rivet_weights_cap),
      std::cref(output_thermodynamics_type),
      std::cref(output_thermodynamics_quantites),
      std::cref(output_thermodynamics_position),
      std::cref(output_thermodynamics_smearing),
      std::cref(output_thermodynamics_onlyParticipants),
      std::cref(lattice_sizes),
      std::cref(lattice_cellNumber),
      std::cref(lattice_origin),
      std::cref(lattice_periodic),
      std::cref(lattice_potentialsAffectThreshold),
      std::cref(potentials_skyrme_skyrmeA),
      std::cref(potentials_skyrme_skyrmeB),
      std::cref(potentials_skyrme_skyrmeTau),
      std::cref(potentials_symmetry_sPot),
      std::cref(potentials_symmetry_gamma),
      std::cref(potentials_vdf_satRhoB),
      std::cref(potentials_vdf_coeffs),
      std::cref(potentials_vdf_powers),
      std::cref(potentials_coulomb_rCut),
      std::cref(forcedThermalization_cellNumber),
      std::cref(forcedThermalization_criticalEDensity),
      std::cref(forcedThermalization_startTime),
      std::cref(forcedThermalization_timestep),
      std::cref(forcedThermalization_algorithm),
      std::cref(forcedThermalization_microcanonical),
      std::cref(forcedThermalization_latticeSizes)};
};

/*!\Userguide
* \page minimum_nonempty_ensembles_
* <hr>
* ### Examples
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
 * ### Example: Configuring the Logging Area
 *
 * To activate different logging levels for different logging areas, change
 the
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
 * This will set all levels to `WARN` verbosity, still asking for
 informational
 * messages of `Main` and `Experiment` areas. Furthermore, `Pythia` debug
 * messages are requested, while any floating point exception message is
 turned
 * off.
 */

/*!\Userguide
 * \page input_collision_term_string_parameters_
 * <hr>
 * ### Example of string paramters customization
 *
 *\verbatim
 Collision_Term:
     Strings: True
     String_Parameters:
         String_Tension: 1.0
         Gluon_Beta: 0.5
         Gluon_Pmin: 0.001
         Quark_Alpha: 2.0
         Quark_Beta: 7.0
         Strange_Supp: 0.16
         Diquark_Supp: 0.036
         Sigma_Perp: 0.42
         StringZ_A_Leading: 0.2
         StringZ_B_Leading: 2.0
         StringZ_A: 2.0
         StringZ_B: 0.55
         String_Sigma_T: 0.5
         Prob_proton_to_d_uu: 0.33
         Separate_Fragment_Baryon: True
         Popcorn_Rate: 0.15
 \endverbatim
 */

/*!\Userguide
 * \page input_collision_term_dileptons_
 * <hr>
 * ### Example of dileptons configuration
 *
 * The following example configures the dilepton production for dileptons
 * originating from resonance decays. In addition, the extended OSCAR2013
 * dilepton output is enabled.
 *
 *\verbatim
 Output:
     Dileptons:
         Format: ["Oscar2013"]
         Extended: True
 Collision_Term:
     Dileptons:
         Decays: True
 \endverbatim
 *
 * <hr>
 * ## Dilepton production in SMASH
 *
 * The treatment of Dilepton Decays is special:
 * - Dileptons are treated via the time integration method, also called
 *   *shining*, as e.g. described in \iref{Schmidt:2008hm}, chapter 2D.
 *   This means that, because dilepton decays are so rare, possible decays are
 *   written in the output at every hadron propagation without ever performing
 *   them. The are weighted with a "shining weight" to compensate for the
 *   over-production.
 * - The shining weight can be found in the weight element of the output.
 * - The shining method is implemented in the DecayActionsFinderDilepton,
 *   which is automatically enabled together with the dilepton output.
 *
 * \anchor input_collision_term_dileptons_note_ \note
 * If you want dilepton decays, you have to modify the *decaymodes.txt* file
 * of your choice, which you then specify as the input with the `-d` command
 * line option. <b>Without this decay modes modification the dilepton output
 * will be empty</b>. Dilepton decays are commented out by default. Therefore,
 * you need to uncomment them. For the N(1520) Dalitz decay, two treatments
 are
 * available: Either by proxy of the \f$\rho N\f$ decay, which is enabled by
 * default (and leads to a dilepton Dalitz decay, if \f$\rho \rightarrow
 * e^+e^-\f$ is also enabled) or as a direct Dalitz decay to \f$e^+e^- N\f$.
 * If using the latter comment-out the \f$\rho N\f$ decay to avoid double
 * counting. The form factor in the direct case, is constant and fixed at the
 * real photon point. Furthermore note, that for dilepton decays, new decay
 * channels can \b not simply be added to the *decaymodes.txt* file. You also
 * have to modify the decay width formulas \c TwoBodyDecayDilepton::width and
 * \c ThreeBodyDecayDilepton::diff_width in *decaytype.cc* file.
 *
 */

/*!\Userguide
 * \page input_collision_term_photons_
 * <hr>
 * ### Example of photons configuration
 *
 * The following example configures the photon production in both binary
 * scatterings and bremsstrahlung processes, where 1000 fractional photons are
 * sampled per single perturbatively produced photon. In addition, the binary
 * photon output is enabled.
 *
 *\verbatim
 Output:
     Photons:
         Format: ["Binary"]
 Collision_Term:
     Photons:
         Fractional_Photons: 1000
         2to2_Scatterings: True
         Bremsstrahlung: True
 \endverbatim
 *
 * <hr>
 * ## Photon production in SMASH
 *
 * Photons are treated perturbatively and are produced from binary
 * scattering processes. Their production follows the framework from Turbide
 * et al. described in \iref{Turbide:2006zz}. Following the perturbative
 * treatment, the produced photons do not contribute to the evolution of the
 * hadronic system. They are rather direcly printed to the photon output.
 * The mechanism for photon production is the following:
 * -# Look for hadronic interactions of particles that are also incoming
 *    particles of a photon process. Currently, the latter include binar
 *    scatterings of \f$ \pi \f$ and \f$ \rho \f$ mesons in the case of
 *    photons from 2-to-2-scatterings or \f$ \pi \f$ scatterings in the case
 *    of bremsstrahlung photons.
 * -# Perform the photon action and write the results to the photon output.
 *    The final state particles are not of interest anymore as they are not
 *    propagated further in the evolution. To account for the probability that
 *    photon processes are significantly less likely than hadronic processes,
 *    the produced photons are weighted according to the ratio of the photon
 *    cross section to the hadronic cross section used to find the
 interaction,
 *    \f[W = \frac{\sigma_\gamma}{\sigma_\mathrm{hadronic}}\;.\f]
 *    This weight can be found in the weight element of the photon output,
 *    denoted as `photon_weight` there.
 * -# Perform the original hadronic action based on which the photon action
 *    was found. Propagate all final states particles throughout the hadronic
 *    evolution as if no photon action had occured.
 *
 * As photons are produced very rarely, a lot of statistics is necessery to
 * yield useful results. Alternatively, it it possible to use fractional
 * photons (see \ref output_content_specific_options_
 * "Content-specific output options" on how to activate them).
 * This means that for each produced photon, \f$ N_{\text{Frac}} \f$
 * photons are actually sampled with different kinematic properties so that
 * more phase space is covered. In case fractional photons are used, the
 * weight for 2-to-2-scatterings is redefined as
 * \f[ W = \frac{\frac{\mathrm{d}\sigma_\gamma}{\mathrm{d}t} \ (t_2 - t_1)}{
 *                     N_\mathrm{frac} \ \sigma_{\mathrm{had}}}. \f]
 *
 * Unlike for binary scatterings, the final state kinematics of bremsstrahlung
 * processes are not entirly defined from the incoming particles. Moreover,
 * the final state momentum of the photon and as well as the scattering angle
 * with respect to the incoming pion collision axis are free parameters whose
 * distribution is encapsulated in the the differential cross section
 * \f$ \frac{\mathrm{d}^2\sigma_\gamma}{\mathrm{d}k\ \mathrm{d} \theta}\f$.
 * For numerical reasons and as the differential cross section can be
 * approximately factorized over the common \f$ k \f$ and
 * \f$ \theta \f$ range, \f$ \frac{\mathrm{d}\sigma_\gamma}{\mathrm{d}k}\f$
 * and \f$ \frac{\mathrm{d}\sigma_\gamma}{\mathrm{d} \theta}\f$ are considered
 * separately. Consequently, the weighting factor in the case of
 bremsstrahlung
 * photons is redefined as:
 * \f[
 * W = \frac{
 *       \sqrt{\frac{\mathrm{d}\sigma_\gamma}{\mathrm{d}k} \ \Delta k \
 *             \frac{\mathrm{d}\sigma_\gamma}{\mathrm{d}\theta}\ \Delta\theta}
 *     }{N_\mathrm{frac}\ \sigma_{\mathrm{had}}}\;,
 * \f]
 * where \f$ \Delta k \f$ and \f$ \Delta\theta \f$ correspond to the
 * available \f$ k \f$ and \f$ \theta \f$ ranges.
 *
 * \note As photons are treated perturbatively, the produced photons are only
 * written to the photon output, but neither to the usual collision output,
 * nor to the particle lists.
 */

/*!\Userguide
 * \page input_modi_collider_
 * <hr>
 * ### Example of heavy-ion collision configuration
 *
 * The following example configures a Cu63-Cu63 collision at
 * \f$\sqrt{s_{NN}}=3.0\,\mathrm{GeV}\f$ with zero impact parameter and Fermi
 * motion taken into consideration. The calculation frame is the default,
 center
 * of velocity, and the nuclei are not deformed. Refer to \ref
 * input_modi_collider_projectile_and_target_ for information about the
 * `Particles` and `Target` sections.
 *
 *\verbatim
 Modi:
     Collider:
         Projectile:
             Particles: {2212: 29, 2112 :34}
         Target:
             Particles: {2212: 29, 2112 :34}
         Sqrtsnn: 3.0
 \endverbatim
 *
 * To further use Fermi motion and allow the first collisions within the
 * projectile or target nucleus, the corresponding options need to be
 activated
 * by means of:
 *\verbatim
         Fermi_Motion: "on"
         Collisions_Within_Nucleus: True
 \endverbatim
 *
 * Additionally, the impact parameter may be specified manually. See
 * \ref input_modi_collider_impact_parameter_ for an example.
 * <hr>
 *
 * \note
 * By default, executing SMASH from the codebase build folder without further
 * specifying the configuration, particles and decay modes files, a collider
 * simulation is set up according to the default _config.yaml_,
 _particles.txt_
 * and _decaymodes.txt_ files located in the _**input**_ directory at the
 * top-level of the codebase. However, changing the _**input**_ directory
 * content will not affect the default SMASH run, unless a clean build folder
 is
 * created over again. This is because the triplet of input files are
 * transformed into another triplet of files into the build directory when
 * `cmake` is run. Hence prefer to use `smash` command line options in case
 you
 * want to refer to possibly modified configuration, particles and decay modes
 * files.\n
 * To run SMASH in the (default) collider setup, execute
 * \verbatim
    ./smash
 \endverbatim
 * from the codebase build folder.
 */

/*!\Userguide
 * \page input_modi_collider_projectile_and_target_
 * <hr>
 * \anchor input_modi_collider_projectile_and_target_ex1_
 * ### p-Pb collisions at the LHC
 *
 * Note that SMASH performs its calculation in the centre-of-velocity and the
 * particles are returned in the centre-of-mass frame. The particles therefore
 * need to be boosted by the rapidity of the centre-of-mass (-0.465 for p-Pb
 * at 5.02TeV).
 * \verbatim
 Modi:
     Collider:
         Calculation_Frame: center of velocity
         Impact:
             Random_Reaction_Plane: True
             Range: [0, 8.5]
         Projectile:
             E_Tot: 1580
             Particles:
                 2212: 82
                 2112: 126
         Target:
             E_Tot: 4000
             Particles:
                 2212: 1
                 2112: 0
 \endverbatim
 *
 * <hr>
 * \anchor input_modi_collider_projectile_and_target_ex2_
 * ### Configuring custom nuclei from external file
 *
 * The following example illustrates how to configure a center-of-mass
 heavy-ion
 * collision with nuclei generated from an external file. The nucleon
 positions
 * are not sampled by SMASH but read in from an external file. The given path
 * and name of the external file are made up and should be defined by the user
 * according to the used file.
 *\verbatim
 Modi:
     Collider:
         Projectile:
             Particles: {2212: 79, 2112: 118}
             Custom:
                 File_Directory: "/home/username/custom_lists"
                 File_Name: "Au197_custom.txt"
         Target:
             Particles: {2212: 79, 2112: 118}
             Custom:
                 File_Directory: "/home/username/custom_lists"
                 File_Name: "Au197_custom.txt"
         Sqrtsnn: 7.7
 \endverbatim
 *
 * The _Au197_custom.txt_ file should be formatted as follows:
 *
 * <div class="fragment">
 * <div class="line"><span class="preprocessor"> 0.20100624   0.11402423
 * -2.40964466   0   0</span></div>
 * <div class="line"><span class="preprocessor"> 1.69072087  -3.21471918
 *  1.06050693   0   1</span></div>
 * <div class="line"><span class="preprocessor">-1.95791109  -3.51483782
 *  2.47294656   1   1</span></div>
 * <div class="line"><span class="preprocessor"> 0.43554894   4.35250733
 *  0.13331011   1   0</span></div>
 * <div class="line"><span class="preprocessor"> ...</span></div>
 * </div>
 *
 * It contains 5 columns (x, y, z, s, c). The first three columns specify the
 * spatial cordinates <b>in fm</b>. The fourth column denotes the spin
 * projection. The fifth contains the charge with 1 and 0 for protons and
 * neutrons respectively. In the example given the first line defines a
 neutron
 * and the second one a proton. Please make sure that your file contains as
 many
 * particles as you specified in the configuration. For the example considered
 * here, the file needs to contain 79 protons and 118 neutrons in the first
 197
 * lines. And the same number in the following 197 lines. The read in nuclei
 are
 * randomly rotated and recentered. Therefore you can run SMASH even if your
 * file does not contain enough nuclei for the number of events you want to
 * simulate as the missing nuclei are generated by rotation of the given
 * configurations.
 *
 * \note
 * SMASH is shipped with an example configuration file to set up a collision
 * with externally generated nucleon positions. This requires a particle list
 to
 * be read in. Both, the configuration file and the particle list, are located
 * in the _**input/custom_nucleus**_ folder at the top-level of SMASH
 codebase.
 * To run SMASH with the provided example configuration and particle list,
 execute
 * \verbatim
    ./smash -i INPUT_DIR/custom_nucleus/config.yaml
 \endverbatim
 * where `INPUT_DIR` needs to be replaced by the path to the input directory
 * at the top-level of SMASH codebase.
 *
 * <hr>
 * \anchor input_modi_collider_projectile_and_target_ex2_
 * ### Configuring a deformed nucleus
 *
 * To configure a fixed target heavy-ion collision with deformed nuclei, whose
 * spherical deformation is explicitly declared, it can be done according to
 the
 * following example. For explanatory (and not physics) reasons, the
 * projectile's Woods-Saxon distribution is initialized automatically and
 * its spherical deformation manually, while the target nucleus is configured
 * just the opposite.
 *\verbatim
 Modi:
     Collider:
         Projectile:
             Particles: {2212: 29, 2112: 34}
             Deformed:
                 # Manually set deformation parameters
                 Automatic: false
                 Beta_2: 0.1
                 Beta_3: 0.2
                 Beta_4: 0.3
                 Orientation:
                     Theta: 0.8
                     Phi: 0.02
                     Psi: 0.13
         Target:
             Particles: {2212: 29, 2112: 34}
             # manually set Woods-Saxon parameters
             Saturation_Density: 0.1968
             Diffusiveness: 0.8
             Radius: 2.0
             Deformed:
                 # Automatically set deformation parameters
                 Automatic: true
                 Orientation:
                     # Randomly rotate nucleus
                     Random_Rotation: true
         E_kin: 1.2
         Calculation_Frame: "fixed target"
 \endverbatim
 */

/*!\Userguide
 * \page input_modi_collider_impact_parameter_
 * <hr>
 * ### Configuring the Impact Parameter
 *
 * The impact parameter can be configured to have a fixed value in the
 * `Collider` subsection of `Modi`. In addition, the initial distance of the
 * nuclei in \f$z\f$-direction is assigned a specific value. This does not
 * affect the time at which the nuclei will collide, but only changes the
 start
 * time of the simulation as the nuclei are further apart when the simulation
 * begins.
 *\verbatim
 Modi:
     Collider:
         Impact:
             Value: 0.1
 \endverbatim
 * The impact parameter may further be sampled within a certain impact
 parameter
 * range. By default, a quadratic distribution is used for the sampling.
 * However, this may be set to `"uniform"` if necessary.
 *\verbatim
 Modi:
     Collider:
         Impact:
             Sample: "quadratic"
             Range: [3.0, 6.0]
 \endverbatim
 * A custom impact parameter distribution based on a set of `Values` and
 * `Yields`, can be configured as follows:
 *\verbatim
 Modi:
     Collider:
         Impact:
             Sample: "custom"
             Values: [0.0, 3.0, 6.0, 9.0]
             Yields: [0.000000, 2.999525, 5.959843, 6.995699]
 \endverbatim
 */

/*!\Userguide
 * \page input_modi_sphere_
 * <hr>
 * ### Configuring a sphere simulation
 *
 * The following example configures an expanding sphere with a radius of 5 fm
 * at a temperature of 200 MeV. The particles are initialized with thermal
 * momenta at a start time of 0 fm. The particle numbers at initialization are
 * 100 \f$ \pi^+ \f$, 100 \f$ \pi^0 \f$, 100 \f$ \pi^- \f$, 50 protons and 50
 * neutrons.
 *
 *\verbatim
 Modi:
     Sphere:
         Radius: 5.0
         Temperature: 0.2
         Initial_Condition: "thermal momenta"
         Start_Time: 0.0
         Init_Multiplicities:
             211: 100
             111: 100
             -211: 100
             2212: 50
             2112: 50
 \endverbatim
 *
 * It is also possible to initialize a sphere based on thermal multiplicities.
 * This is done via
 *\verbatim
 Modi:
     Sphere:
         Radius: 10.0
         Temperature: 0.2
         Use_Thermal_Multiplicities: True
 \endverbatim
 *
 * If one wants to simulate a jet in the hadronic medium, this can be done by
 * using the following configuration setup:
 *\verbatim
 Modi:
     Sphere:
         Radius: 10.0
         Temperature: 0.2
         Use_Thermal_Multiplicities: True
         Jet:
             Jet_PDG: 211
             Jet_Momentum: 100.0
\endverbatim
 *
 * \note
 * SMASH is shipped with an example configuration file to set up an expanding
 * sphere simulation initialized with predefined initial particle
 * multiplicities. This file is located in the _**input/sphere**_ folder at
 * the top-level of SMASH codebase. To run SMASH with the provided example
 * configuration for the sphere system, execute
 * \verbatim
    ./smash -i INPUT_DIR/sphere/config.yaml
 \endverbatim
 * where `INPUT_DIR` needs to be replaced by the path to the input directory
 * at the top-level of SMASH codebase.
 *
 */

/*!\Userguide
 * \page input_modi_box_
 * <hr>
 * ### Configuring a Box Simulation
 *
 * The following example configures an infinite matter simulation in a Box with
 * 10 fm cube length at a temperature of 200 MeV. The particles are initialized
 * with thermal momenta at a start time of 10 fm. The particle numbers at
 * initialization are 100 \f$ \pi^+ \f$, 100 \f$ \pi^0 \f$, 100 \f$ \pi^- \f$,
 * 50 protons and 50 neutrons.
 *
 *\verbatim
 Modi:
     Box:
         Length: 10.0
         Temperature: 0.2
         Initial_Condition: "thermal momenta"
         Start_Time: 10.0
         Init_Multiplicities:
             211: 100
             111: 100
             -211: 100
             2212: 50
             2112: 50
 \endverbatim
 * On the contrary, it is also possible to initialize a thermal box based on
 * thermal multiplicities. This is done via
 *\verbatim
 Modi:
     Box:
         Length: 10.0
         Temperature: 0.2
         Use_Thermal_Multiplicities: True
         Initial_Condition: "thermal momenta"
         Baryon_Chemical_Potential: 0.0
         Strange_Chemical_Potential: 0.0
         Charge_Chemical_Potential: 0.0
         Account_Resonance_Widths: True
 \endverbatim
 *
 * If one wants to simulate a jet in the hadronic medium, this can be done
 * by using the following configuration setup:
 *\verbatim
 Modi:
     Box:
         Length: 10.0
         Temperature: 0.2
         Use_Thermal_Multiplicities: True
         Initial_Condition: "thermal momenta"
         Jet:
             Jet_PDG: 211
             Jet_Momentum: 100.0
\endverbatim
 *
 * \note
 * The box modus is most useful for infinite matter simulations with thermal and
 * chemical equilibration and detailed balance. Detailed balance can however not
 * be conserved if 3-body decays (or higher) are performed. To yield useful
 * results applying a SMASH box simulation, it is therefore necessary to modify
 * the provided default _particles.txt_ and _decaymodes.txt_ files by removing
 * 3-body and higher order decays from the decay modes file and all
 * corresponding particles that can no longer be produced from the particles
 * file. In addtion, strings need to be turned off, since they also break
 * detailed balance due to lacking backreactions.\n\n
 * SMASH is shipped with example files (_config.yaml_, _particles.txt_ and
 * _decaymodes.txt_) meeting the above mentioned requirements to set up an
 * infinite matter simulation. These files are located in the _**input/box**_
 * folder at the top-level of SMASH codebase. To run SMASH with the provided
 * example configuration for the box system, execute
 * \n
 * \verbatim
    ./smash -i INPUT_DIR/box/config.yaml\
            -p INPUT_DIR/box/particles.txt\
            -d INPUT_DIR/box/decaymodes.txt
 \endverbatim
 * where `INPUT_DIR` needs to be replaced by the path to the input directory
 * at the top-level of SMASH codebase.
 */

/*!\Userguide
 * \page input_modi_list_
 * <hr>
 * ### Configuring an afterburner simulation
 *
 * The following example sets up an afterburner simulation for a set of particle
 * files located in _**particle_lists_in**_ folder. The files are named as
 * _event10_, _event11_, etc. (the first being number 10 is specified by the key
 * `Shift_Id`). SMASH is run once for each event in the folder.
 * \verbatim
 Modi:
     List:
         File_Directory: "particle_lists_in"
         File_Prefix: "event"
         Shift_Id: 10
 \endverbatim
 *
 * <hr>
 * ## Some information about the structure of input particle file
 *
 * This is how an input particle file might look like:
 * <div class="fragment">
 * <div class="line"><span class="preprocessor">#!OSCAR2013 particle_lists
 * t x y z mass p0 px py pz pdg ID charge</span></div>
 * <div class="line"><span class="preprocessor">\# Units: fm fm fm fm
 * GeV GeV GeV GeV GeV none none none</span></div>
 * <div class="line"><span class="preprocessor">0.1 6.42036 1.66473 9.38499
 * 0.138 0.232871 0.116953 -0.115553 0.090303 111 0 0</span></div>
 * <div class="line"><span class="preprocessor">\# event 0 end</span></div>
 * <div class="line"><span class="preprocessor">\# event 1</span></div>
 * <div class="line"><span class="preprocessor">0.1 6.42036 1.66473 9.38499
 * 0.138 0.232871 0.116953 -0.115553 0.090303 111 0 0</span></div>
 * <div class="line"><span class="preprocessor">\# event 1 end</span></div>
 * </div>
 * Each colum contains the described quantities. In particular, in the example
 * above, one \f$\pi^0\f$ with spatial coordinates
 * \f[(t, x, y, z) = (0.1, 6.42036, 1.66473, 9.38499)\,\mathrm{fm}\f]
 * and and 4-momenta
 * \f[(p_0,p_x,p_y,p_z)=(0.232871,0.116953,-0.115553,0.090303)\,\mathrm{GeV}\f]
 * with mass = 0.138 GeV, pdg = 111, id = 0 and charge 0 will be initialized for
 * the first event (and also for the second event).
 *
 * \note
 * SMASH is shipped with an example configuration file to set up an afterburner
 * simulation by means of the list modus. This also requires a particle list to
 * be read in. Both, the configuration file and the particle list, are located
 * in the _**input/list**_ folder at the top-level of SMASH codebase. To run
 * SMASH with the provided example configuration and particle list, execute
 * \verbatim
    ./smash -i INPUT_DIR/list/config.yaml
 \endverbatim
 * where `INPUT_DIR` needs to be replaced by the path to the input directory
 * at the top-level of SMASH codebase.
 */

/*!\Userguide
 * \page input_lattice_
 * <hr>
 * ### Configuring the Lattice
 *
 * The following example configures the lattice with the origin in (0,0,0), 20
 * cells of 10 fm size in each direction and with periodic boundary conditions.
 * The potential effects on the thresholds are taken into consideration. Note
 * that, as the origin is by definition the left down near corner of the cell,
 * center is located at (5, 5, 5).
 *\verbatim
 Lattice:
     Origin: [0.0, 0.0, 0.0]
     Sizes: [10.0, 10.0, 10.0]
     Cell_Number: [20, 20, 20]
     Periodic: True
     Potentials_Affect_Thresholds: True
 \endverbatim
 * In case of `"Collider"`, `"Box"`, and `"Sphere"` <tt>\ref key_gen_modus_
 * "Modus"</tt> there is also an option to set up lattice automatically. For
 * example, for `"Collider"` modus
 *\verbatim
 Lattice:
 \endverbatim
 * sets up a lattice that (heuristically) covers causally allowed particle
 * positions until the <tt>\ref key_gen_end_time_ "End_Time"</tt> of the
 * simulation. The lattice may also be automatically contracted in z direction
 * depending on the chosen way of density calculation.
 *
 * Another example for `"Box"` modus:
 *\verbatim
 Lattice:
     Cell_Number: [20, 20, 20]
 \endverbatim
 * This sets up a periodic lattice that matches box sizes.
 */

/*!\Userguide
 * \page input_forced_thermalization_
 * <hr>
 * ### Configuring forced thermalization
 *
 * The following example activates forced thermalization in cells in which the
 * energy density is above 0.3 GeV/fm³. The lattice is initialized with 21
 * cells in x and y direction and 101 cells in z-direction. The lattice size is
 * 20 fm in x and y direction and 50 fm in z-direction. The thermalization is
 * applied only for times later than 10 fm with a timestep of 1 fm/c. The
 * sampling is done according to the "biased BF" algorithm.
 *\verbatim
 Forced_Thermalization:
     Lattice_Sizes: [20.0, 20.0, 50.0]
     Cell_Number: [21, 21, 101]
     Critical_Edens: 0.3
     Start_Time: 10.0
     Timestep: 1.0
     Algorithm: "biased BF"
 \endverbatim
 */

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_VALIDATION_H_