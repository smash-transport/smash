/*
 *
 *    Copyright (c) 2014-2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_CONFIGURATION_H_
#define SRC_INCLUDE_SMASH_CONFIGURATION_H_

#include <array>
#include <filesystem>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "yaml-cpp/yaml.h"

#include "density.h"
#include "forwarddeclarations.h"

namespace YAML {

/**
 * Convert from YAML::Node to SMASH-readable (C++) format and vice versa.
 *
 * \tparam T Type of the values (could be any data type that
 * needs conversion).
 */
template <typename T>
struct convert {
  /**
   * Serialization: Converts x (of any type) to a YAML::Node. To do this,
   * the type of x needs first be cast to a string.
   *
   * \param[in] x Value that is to be converted to a YAML::Node.
   * \return YAML node
   */
  static Node encode(const T &x) { return Node{static_cast<std::string>(x)}; }

  /**
   * Deserialization: Converts a YAML::Node to any SMASH-readable data type and
   * returns whether or not this operation was successful.
   *
   * \param[in] node YAML::Node that is to be converted.
   * \param[in] x Value that the YAML:Node is cast to.
   * \return True in case conversion was successful.
   */
  static bool decode(const Node &node, T &x) {
    if (!node.IsScalar()) {
      return false;
    } else {
      x = static_cast<T>(node.Scalar());
      return true;
    }
  }
};
}  // namespace YAML

namespace smash {

/*!\Userguide
 * \page inputparticles Particles
 *
 * The particles available to SMASH are defined in `particles.txt`, which is
 * located in '$SMASH_SRC_DIRECTORY/input'. If you want to modify and use this
 * file to set up SMASH, execute
 * ```
 * ./smash -p $SMASH_SRC_DIRECTORY/input/particles.txt
 * ```
 * in the '$SMASH_SRC_DIRECTORY/build' directory. \n
 *
 * The particles are
 * given as a table with the particles properties in different columns. Note,
 * that these columns may be separated by an arbitrary number of spaces:
 * ```
 * <name> <mass in GeV> <width in GeV> <parity> <PDG codes>
 * ```
 *
 * The name has to be a unique UTF-8 string. Conventionally, unicode names are
 * used in SMASH to make the file more readable and generate prettier output. It
 * is possible to only specify the isospin multiplet and SMASH will fill in the
 * properties of the components of the multiplet assuming isospin symmetry. The
 * names generated this way will have the charges appended to the multiplet name
 * using the unicode characters `⁻`, `⁰` and `⁺`. This is appropriate for almost
 * all particles. Anti particles do not have to be specified explicitely.
 *
 * The pole mass and the on-shell width of the particle or multiplet have to be
 * specified as floating point numbers in GeV.
 *
 * The parity has to be either `+` or `-`.
 *
 * The PDG codes are following the [numbering
 * scheme](http://pdg.lbl.gov/2018/mcdata/mc_particle_id_contents.html)
 * specified by the PDG, which depends on the quantum numbers of the particles.
 * For SMASH, it is important that the quark content in the PDG code is
 * correctly specified. Other than that, deviations from the numbering scheme
 * have no effect in SMASH. If the name represents a multiplet, there has to be
 * a PDG code for all multiplet members, except for anti particles.
 *
 * For example, to define all three pions (π⁻, π⁰, π⁺), it is sufficient to
 * specify the π multiplet using the following line in `particles.txt`, where
 * the 4th column contains the PDG number of the neutral and the 5th PDG number
 * of the charged state:
 * ```
 * π  0.138  7.7e-9  111  211
 * ```
 *
 * It is also possible to only specify a specific member of the multiplet. In
 * this case, the charge has to be given as a suffix in the name using the UTF-8
 * unicode characters `⁻`, `⁰` and `⁺`. For example, the properties of the
 * electron can be specified like this:
 * ```
 * e⁻  0.000511  0  11
 * ```
 *
 * Comments can be added to `particles.txt` using the `#` character. Everything
 * after `#` until the end of the line is ignored.
 *
 * Note that some reactions in SMASH are parametrized and require specific
 * particles in the final state. When such a reaction happens and the required
 * particle is not defined, SMASH will crash.
 *
 * If you specify an incorrect value, SMASH will print an error similar to the
 * following:
 * ```
 * Failed to convert the input string to the expected data types.
 * ```
 *
 * When running a box simulation in which detailed balance is expected to be
 * conserved, the particles file will need to be modified. See \ref
 * input_modi_box_ for further information.
 */

/*!\Userguide
 * \page inputdecaymodes Decay Modes
 *
 *All possible decays and resonance formations in SMASH are provided by the
 * `decaymodes.txt` file, which is
 * located in '$SMASH_SRC_DIRECTORY/input'. If you want to modify and use this
 * file to set up SMASH, execute
 * ```
 * ./smash -d $SMASH_SRC_DIRECTORY/input/decaymodes.txt
 * ```
 * in the '$SMASH_SRC_DIRECTORY/build' directory. \n
 *
 * The decaymodes are formatted in blocks of the following format:
 * ```
 * <name of decaying particle>
 * <branching ratio> <angular momentum L> <names of decay products>
 * <branching ratio> <angular momentum L> <names of decay products>
 * ...
 * ```
 * The blocks have to be separated by at least one empty line.
 *
 * The names have to be the ones defined in `particles.txt` (see \ref
 * inputparticles). If multiplet names are used, the other branching ratios are
 * generated by SMASH assuming isospin symmetry. Note that currently decay
 * channels can only be specified for whole multiplets; individual particles
 * can however still be used in a decay channel as specific daughters.
 *
 * The branching ratios are given as a floating point number. If the branching
 * ratios in one block do not add up to 1, they are automatically normalized by
 * SMASH.
 *
 * The angular momentum of the decay channel has to be specified as an integer.
 *
 * The names of two or three decay products have to be given for each channel.
 * Note that the SMASH defaults avoid three-body decays, because they break
 * detailed balance due to the lack of 3-to-1 reactions in SMASH.
 *
 * For example, the following lines are enough to specify all possible decays of
 * the N(1440) resonance multiplet:
 * ```
 * N(1440)
 * 0.60   1  N π
 * 0.24   1  Δ π
 * 0.16   0  N σ
 * ```
 * For decays violating isospin symmetry, it is possible to specify the members
 * of the multiplets in the final state explicitely:
 * ```
 * φ
 * 0.489   1  K⁺ K̅⁻
 * 0.342   1  K⁰ K̅⁰
 * ```
 *
 * It is possible to add comments to `decaymodes.txt` using the `#` character.
 * Everything after `#` until the end of the line is ignored.
 *
 * Note that SMASH has an internal width cut-off (currently 10 keV), below which
 * particles cannot decay, even if decays are specified in `decaymodes.txt`.
 * This only changes when Include_Weak_And_EM_Decays_At_The_End is enabled, then
 * all decays are considered in the final decays, no matter their decay width .
 *
 * Note further, that the decaymodes file will need to be modified when running
 * a box simulation in which detailed balance is expected to be conserved. See
 * \ref input_modi_box_ for further information.
 */

/**
 * Interface to the SMASH configuration files.
 *
 * The configuration is created from a YAML file and then stores a nested map of
 * maps (normally a tree, but YAML allows it to be cyclic - even though we don't
 * want that feature).
 *
 * For the typical usage in SMASH one needs to read the value once. In that
 * case, use the Configuration::take function, for example:
 * \code
 * double sigma = config.take({"General", "SIGMA"});
 * \endcode
 * Note the curly braces in the function call. It is a std::initializer_list of
 * strings. This allows an arbitrary nesting depth via the same function.
 * But as a consequence the keys must all be given as constant strings at
 * compile time.
 *
 * If you need to access the configuration values from a run-time string you can
 * use Configuration::operator[]. This returns a Configuration object that
 * references the respective sub-tree.
 *
 * By taking values (instead of just reading), the configuration object should
 * be empty at the end of the initialization. If the object is not empty, SMASH
 * will print a warning (using Configuration::unused_values_report). This can be
 * important for the user to discover typos in his configuration file (or
 * command line parameters).
 */
class Configuration {
 public:
  /**
   * \ingroup exception
   * Thrown when the types in the config file and C++ don't match.
   */
  struct IncorrectTypeInAssignment : public std::runtime_error {
    using std::runtime_error::runtime_error;
  };
  /**
   * \ingroup exception
   * Thrown for YAML parse errors.
   */
  struct ParseError : public std::runtime_error {
    using std::runtime_error::runtime_error;
  };

  /**
   * \ingroup exception
   * Thrown if the file does not exist.
   */
  struct FileDoesNotExist : public std::runtime_error {
    using std::runtime_error::runtime_error;
  };

  /**
   * Return type of Configuration::take that automatically determines the target
   * type.
   *
   * This class is an implementation detail of Configuration and can be ignored
   * by users of Configuration.
   */
  class Value {
    friend class Configuration;

    /// a YAML leaf node \todo(steinberg) What is that?
    const YAML::Node node_;
    /// The key to be interpreted
    const char *const key_;

    /**
     * Construct the Value wrapper from a YAML::Node.
     *
     * \note This constructor must be implicit, otherwise it's impossible to
     * return an rvalue Value object - because the copy constructor is deleted.
     */
    Value(const YAML::Node &n, const char *key) : node_(n), key_(key) {
      if (!(n.IsScalar() || n.IsSequence() || n.IsMap())) {
        std::stringstream err;
        err << "Configuration value for \"" << key
            << "\" is missing or invalid";
        throw std::runtime_error(err.str());
      }
    }

   public:
    /// If you want to copy this you're doing it wrong
    Value(const Value &) = delete;
    /// If you want to copy this you're doing it wrong
    Value &operator=(const Value &) = delete;

    /**
     * Convert the value to the type of the supplied argument.
     *
     * The argument itself is not used other than to determine its type. This
     * function is necessary because in some situations the overload resolution
     * rules lead to the correct conversion becoming hidden. Then you'll see a
     * compiler error with a list of ambiguous constructor calls as candidates.
     * Use this function as a workaround.
     * Example:
     * \code
     * // this doesn't compile:
     * const PdgCode code0(config.take({"key"}));
     * // this compiles (because PdgCode::operator= is not overloaded), but note
     * // that it cannot be used in constructor initializer lists:
     * const PdgCode code1 = config.take({"key"});
     *
     * // Thus, for class member variables use the following pattern:
     * class X {
     *  public:
     *   X() : code_(config.take({"key"}).convert_for(code_)) {}
     *
     *  private:
     *   const PdgCode code_;
     * };
     * \endcode
     */
    template <typename T>
    T convert_for(const T &) const {
      return *this;
    }

    /**
     * This function determines the type it is assigned to and calls
     * YAML::Node::as<T>() with this type.
     *
     * This makes reading values more convenient than calling as<type>()
     * explicitly.
     * \throw IncorrectTypeInAssignment
     */
    template <typename T>
    operator T() const {
      try {
        return node_.as<T>();
      } catch (YAML::TypedBadConversion<T> &e) {
        throw IncorrectTypeInAssignment(
            "The value for key \"" + std::string(key_) +
            "\" cannot be converted to the requested type.");
      }
    }

    /**
     * Check conversion exceptions.
     *
     * \throw IncorrectTypeInAssignment in case type conversion failed.
     */
    template <typename T>
    operator std::vector<T>() const {
      try {
        return node_.as<std::vector<T>>();
      } catch (YAML::TypedBadConversion<T> &e) {
        throw IncorrectTypeInAssignment(
            "One of the values in the sequence for key \"" + std::string(key_) +
            "\" failed to convert to the requested type. E.g. [1 2] is a "
            "sequence of one string \"1 2\" and [1, 2] is a sequence of two "
            "integers. Often there is just a comma missing in the config "
            "file.");
      } catch (YAML::TypedBadConversion<std::vector<T>> &e) {
        throw IncorrectTypeInAssignment(
            "The value for key \"" + std::string(key_) +
            "\" cannot be converted to the requested type. A sequence was "
            "expected but apparently not found.");
      }
    }

    /**
     * Cast array of keys to a std::array of length N.
     *
     * \return Array of std::array type.
     * \throw IncorrectTypeInAssignment in case the number of keys does not
     * match the length of the newly generated array.
     */
    template <typename T, size_t N>
    operator std::array<T, N>() const {
      const std::vector<T> vec = operator std::vector<T>();
      const size_t n_read = vec.size();
      // Alert if size does not match
      if (n_read != N) {
        throw IncorrectTypeInAssignment("Wrong number of values in array \"" +
                                        std::string(key_) + "\". Expected " +
                                        std::to_string(N) +
                                        " values,"
                                        " found " +
                                        std::to_string(n_read) + ".");
      }
      std::array<T, N> arr;
      std::copy_n(vec.begin(), N, arr.begin());
      return arr;
    }

    /**
     * Set ReactionBitSet from configuration values.
     *
     * \return ReactionBitSet with all included reaction types.
     * \throw IncorrectTypeInAssignment in case a reaction type that is not
     * available is provided as a configuration value.
     */
    operator ReactionsBitSet() const {
      const std::vector<std::string> v = operator std::vector<std::string>();
      ReactionsBitSet s;
      for (const auto &x : v) {
        if (x == "All") {
          s.set();
          break;
        } else if (x == "Elastic") {
          s.set(IncludedReactions::Elastic);
        } else if (x == "NN_to_NR") {
          s.set(IncludedReactions::NN_to_NR);
        } else if (x == "NN_to_DR") {
          s.set(IncludedReactions::NN_to_DR);
        } else if (x == "KN_to_KN") {
          s.set(IncludedReactions::KN_to_KN);
        } else if (x == "KN_to_KDelta") {
          s.set(IncludedReactions::KN_to_KDelta);
        } else if (x == "Strangeness_exchange") {
          s.set(IncludedReactions::Strangeness_exchange);
        } else if (x == "NNbar") {
          s.set(IncludedReactions::NNbar);
        } else if (x == "PiDeuteron_to_NN") {
          s.set(IncludedReactions::PiDeuteron_to_NN);
        } else if (x == "PiDeuteron_to_pidprime") {
          s.set(IncludedReactions::PiDeuteron_to_pidprime);
        } else if (x == "NDeuteron_to_Ndprime") {
          s.set(IncludedReactions::NDeuteron_to_Ndprime);
        } else {
          throw IncorrectTypeInAssignment(
              "The value for key \"" + std::string(key_) +
              "\" should be \"All\", \"Elastic\", \"NN_to_NR\", \"NN_to_DR\","
              "\"KN_to_KN\", \"KN_to_KDelta\", \"PiDeuteron_to_NN\", "
              "\"PiDeuteron_to_pidprime\", \"NDeuteron_to_Ndprime\", "
              "\"Strangeness_exchange\" or "
              "\"NNbar\", or any combination of these.");
        }
      }
      return s;
    }

    /**
     * Set MultiParticleReactionsBitSet from configuration values.
     *
     * \return MultiParticleReactionsBitSet with all included reaction types.
     * \throw IncorrectTypeInAssignment in case a reaction type that is not
     * available is provided as a configuration value.
     */
    operator MultiParticleReactionsBitSet() const {
      const std::vector<std::string> v = operator std::vector<std::string>();
      MultiParticleReactionsBitSet s;
      for (const auto &x : v) {
        if (x == "All") {
          s.set();
          break;
        } else if (x == "Meson_3to1") {
          s.set(IncludedMultiParticleReactions::Meson_3to1);
        } else if (x == "Deuteron_3to2") {
          s.set(IncludedMultiParticleReactions::Deuteron_3to2);
        } else if (x == "NNbar_5to2") {
          s.set(IncludedMultiParticleReactions::NNbar_5to2);
        } else if (x == "A3_Nuclei_4to2") {
          s.set(IncludedMultiParticleReactions::A3_Nuclei_4to2);
        } else {
          throw IncorrectTypeInAssignment(
              "The value for key \"" + std::string(key_) +
              "\" should be \"All\", \"Meson_3to1\", "
              "\"Deuteron_3to2\" or \"NNbar_5to2\", "
              "\"A3_Nuclei_4to2\", or any combination of "
              "these.");
        }
      }
      return s;
    }

    /**
     * Set thermodynamic quantity from configuration values.
     *
     * \return Set of thermodynamic quantity.
     * \throw IncorrectTypeInAssignment in case a thermodynamic quantity that is
     * not available is provided as a configuration value.
     */
    operator std::set<ThermodynamicQuantity>() const {
      const std::vector<std::string> v = operator std::vector<std::string>();
      std::set<ThermodynamicQuantity> s;
      for (const auto &x : v) {
        if (x == "rho_eckart") {
          s.insert(ThermodynamicQuantity::EckartDensity);
        } else if (x == "tmn") {
          s.insert(ThermodynamicQuantity::Tmn);
        } else if (x == "tmn_landau") {
          s.insert(ThermodynamicQuantity::TmnLandau);
        } else if (x == "landau_velocity") {
          s.insert(ThermodynamicQuantity::LandauVelocity);
        } else if (x == "j_QBS") {
          s.insert(ThermodynamicQuantity::j_QBS);
        } else {
          throw IncorrectTypeInAssignment(
              "The value for key \"" + std::string(key_) +
              "\" should be \"rho_eckart\", \"tmn\""
              ", \"tmn_landau\", \"landau_velocity\" or \"j_QBS\".");
        }
      }
      return s;
    }

    /**
     * Set calculation frame from configuration values.
     *
     * \return string of calculation frame.
     * \throw IncorrectTypeInAssignment in case a calculation frame that is
     * not available is provided as a configuration value.
     */
    operator CalculationFrame() const {
      const std::string s = operator std::string();
      if (s == "center of velocity") {
        return CalculationFrame::CenterOfVelocity;
      }
      if (s == "center of mass") {
        return CalculationFrame::CenterOfMass;
      }
      if (s == "fixed target") {
        return CalculationFrame::FixedTarget;
      }
      throw IncorrectTypeInAssignment(
          "The value for key \"" + std::string(key_) +
          "\" should be \"center of velocity\" or \"center of mass\" "
          "or \"fixed target\".");
    }

    /**
     * (De-)Activate Fermi motion from configuration values.
     *
     * \return Fermi motion setup.
     * \throw IncorrectTypeInAssignment in case a Fermi motion value that is
     * not available is provided as a configuration value.
     */
    operator FermiMotion() const {
      const std::string s = operator std::string();
      if (s == "off") {
        return FermiMotion::Off;
      }
      if (s == "on") {
        return FermiMotion::On;
      }
      if (s == "frozen") {
        return FermiMotion::Frozen;
      }
      throw IncorrectTypeInAssignment(
          "The value for key \"" + std::string(key_) +
          "\" should be \"off\" or \"on\" or \"frozen\".");
    }

    /**
     * Set density type from configuration values.
     *
     * \return Density type.
     * \throw IncorrectTypeInAssignment in case a density type that is
     * not available is provided as a configuration value.
     */
    operator DensityType() const {
      const std::string s = operator std::string();
      if (s == "hadron") {
        return DensityType::Hadron;
      }
      if (s == "baryon") {
        return DensityType::Baryon;
      }
      if (s == "baryonic isospin") {
        return DensityType::BaryonicIsospin;
      }
      if (s == "pion") {
        return DensityType::Pion;
      }
      if (s == "total isospin") {
        return DensityType::Isospin3_tot;
      }
      if (s == "none") {
        return DensityType::None;
      }
      throw IncorrectTypeInAssignment("The value for key \"" +
                                      std::string(key_) +
                                      "\" should be \"hadron\" or \"baryon\" "
                                      "or \"baryonic isospin\" or \"pion\" "
                                      "or \"none\".");
    }

    /**
     * Set expansion mode from configuration values.
     *
     * \return Expansion mode.
     * \throw IncorrectTypeInAssignment in case an expansion mode that is
     * not available is provided as a configuration value.
     */
    operator ExpansionMode() const {
      const std::string s = operator std::string();
      if (s == "NoExpansion") {
        return ExpansionMode::NoExpansion;
      }
      if (s == "MasslessFRW") {
        return ExpansionMode::MasslessFRW;
      }
      if (s == "MassiveFRW") {
        return ExpansionMode::MassiveFRW;
      }
      if (s == "Exponential") {
        return ExpansionMode::Exponential;
      }
      throw IncorrectTypeInAssignment(
          "The value for key \"" + std::string(key_) +
          "\" should be \"NoExpansion\", \"MasslessFRW\"," +
          "\"MassiveFRW\" or \"Exponential\".");
    }

    /**
     * Set DerivativesMode.
     */
    operator DerivativesMode() const {
      const std::string s = operator std::string();
      if (s == "Covariant Gaussian") {
        return DerivativesMode::CovariantGaussian;
      }
      if (s == "Finite difference") {
        return DerivativesMode::FiniteDifference;
      }
      if (s == "Off") {
        return DerivativesMode::Off;
      }
      throw IncorrectTypeInAssignment(
          "The value for key \"" + std::string(key_) +
          "\" should be \"Covariant Gaussian\", \"Finite difference\"," +
          " or \"Off\".");
    }

    /**
     * Set FieldDerivatives mode.
     */
    operator FieldDerivativesMode() const {
      const std::string s = operator std::string();
      if (s == "Chain Rule") {
        return FieldDerivativesMode::ChainRule;
      }
      if (s == "Direct") {
        return FieldDerivativesMode::Direct;
      }
      throw IncorrectTypeInAssignment(
          "The value for key \"" + std::string(key_) +
          "\" should be \"Chain Rule\" or \"Direct\".");
    }

    /**
     * Set SmearingMode.
     */
    operator SmearingMode() const {
      const std::string s = operator std::string();
      if (s == "Covariant Gaussian") {
        return SmearingMode::CovariantGaussian;
      }
      if (s == "Discrete") {
        return SmearingMode::Discrete;
      }
      if (s == "Triangular") {
        return SmearingMode::Triangular;
      }
      throw IncorrectTypeInAssignment(
          "The value for key \"" + std::string(key_) +
          "\" should be \"Covariant Gaussian\", \"Discrete\"," +
          " or \"Triangular\".");
    }

    /**
     * Set time step mode from configuration values.
     *
     * \return time step mode.
     * \throw IncorrectTypeInAssignment in case a time step mode that is
     * not available is provided as a configuration value.
     */
    operator TimeStepMode() const {
      const std::string s = operator std::string();
      if (s == "None") {
        return TimeStepMode::None;
      }
      if (s == "Fixed") {
        return TimeStepMode::Fixed;
      }
      throw IncorrectTypeInAssignment("The value for key \"" +
                                      std::string(key_) +
                                      "\" should be \"None\" or \"Fixed\".");
    }

    /**
     * Set initial condition for box setup from configuration values.
     *
     * \return Initial condition for box setup.
     * \throw IncorrectTypeInAssignment in case an initial conditions that is
     * not available is provided as a configuration value.
     */
    operator BoxInitialCondition() const {
      const std::string s = operator std::string();
      if (s == "thermal momenta") {
        return BoxInitialCondition::ThermalMomentaBoltzmann;
      }
      if (s == "thermal momenta quantum") {
        return BoxInitialCondition::ThermalMomentaQuantum;
      }
      if (s == "peaked momenta") {
        return BoxInitialCondition::PeakedMomenta;
      }
      throw IncorrectTypeInAssignment(
          "The value for key \"" + std::string(key_) +
          "\" should be \"thermal momenta\", \"thermal momenta quantum\", " +
          "or \"peaked momenta\".");
    }

    /**
     * Set initial condition for sphere setup from configuration values.
     *
     * \return Initial condition for sphere setup.
     * \throw IncorrectTypeInAssignment in case an initial conditions that is
     * not available is provided as a configuration value.
     */
    operator SphereInitialCondition() const {
      const std::string s = operator std::string();
      if (s == "thermal momenta") {
        return SphereInitialCondition::ThermalMomentaBoltzmann;
      }
      if (s == "thermal momenta quantum") {
        return SphereInitialCondition::ThermalMomentaQuantum;
      }
      if (s == "IC_ES") {
        return SphereInitialCondition::IC_ES;
      }
      if (s == "IC_1M") {
        return SphereInitialCondition::IC_1M;
      }
      if (s == "IC_2M") {
        return SphereInitialCondition::IC_2M;
      }
      if (s == "IC_Massive") {
        return SphereInitialCondition::IC_Massive;
      }
      throw IncorrectTypeInAssignment(
          "The value for key \"" + std::string(key_) +
          "\" should be \"thermal momenta\", \"thermal momenta quantum\", " +
          "\"IC_ES\", \"IC_1M\", \"IC_2M\" or" + "\"IC_Massive\".");
    }

    /**
     * Set treatment of N-Nbar reactions from configuration values.
     *
     * \return N-Nbar treatment.
     * \throw IncorrectTypeInAssignment in case an N-Nbar treatment that is
     * not available is provided as a configuration value.
     */
    operator NNbarTreatment() const {
      const std::string s = operator std::string();
      if (s == "no annihilation") {
        return NNbarTreatment::NoAnnihilation;
      }
      if (s == "resonances") {
        return NNbarTreatment::Resonances;
      }
      if (s == "two to five") {
        return NNbarTreatment::TwoToFive;
      }
      if (s == "strings") {
        return NNbarTreatment::Strings;
      }
      throw IncorrectTypeInAssignment(
          "The value for key \"" + std::string(key_) + "\" should be " +
          "\"no annihilation\", \"resonances\", \"two to five\" or " +
          " \"strings\".");
    }

    /**
     * Set cross-section sampling method from configuration values.
     *
     * \return Sampling method of cross-section.
     * \throw IncorrectTypeInAssignment in case a sampling method that is
     * not available is provided as a configuration value.
     */
    operator Sampling() const {
      const std::string s = operator std::string();
      if (s == "quadratic") {
        return Sampling::Quadratic;
      }
      if (s == "custom") {
        return Sampling::Custom;
      }
      if (s == "uniform") {
        return Sampling::Uniform;
      }
      throw IncorrectTypeInAssignment(
          "The value for key \"" + std::string(key_) +
          "\" should be \"quadratic\", \"uniform\" or \"custom\".");
    }

    /**
     * Set algorithm for forced thermalization from configuration values.
     *
     * \return Algorithm for forced thermalization.
     * \throw IncorrectTypeInAssignment in case a thermalization algorithm that
     * is not available is provided as a configuration value.
     */
    operator ThermalizationAlgorithm() const {
      const std::string s = operator std::string();
      if (s == "mode sampling") {
        return ThermalizationAlgorithm::ModeSampling;
      }
      if (s == "biased BF") {
        return ThermalizationAlgorithm::BiasedBF;
      }
      if (s == "unbiased BF") {
        return ThermalizationAlgorithm::UnbiasedBF;
      }
      throw IncorrectTypeInAssignment(
          "The value for key \"" + std::string(key_) +
          "\" should be \"mode sampling\", \"biased BF\" or \"unbiased BF\".");
    }

    /**
     * Set collision criterion from configuration values.
     *
     * \return CollisionCriterion.
     * \throw IncorrectTypeInAssignment in case an collision criterion that is
     * not available is provided as a configuration value.
     */
    operator CollisionCriterion() const {
      const std::string s = operator std::string();
      if (s == "Geometric") {
        return CollisionCriterion::Geometric;
      }
      if (s == "Stochastic") {
        return CollisionCriterion::Stochastic;
      }
      if (s == "Covariant") {
        return CollisionCriterion::Covariant;
      }
      throw IncorrectTypeInAssignment(
          "The value for key \"" + std::string(key_) + "\" should be " +
          "\"Geometric\", \"Stochastic\" " + "or \"Covariant\".");
    }

    /**
     * Set OutputOnlyFinal for particles output from configuration values.
     *
     * \return OutputOnlyFinal.
     * \throw IncorrectTypeInAssignment in case only_final value that is
     * not available is provided as a configuration value.
     */
    operator OutputOnlyFinal() const {
      const std::string s = operator std::string();
      if (s == "Yes") {
        return OutputOnlyFinal::Yes;
      }
      if (s == "No") {
        return OutputOnlyFinal::No;
      }
      if (s == "IfNotEmpty") {
        return OutputOnlyFinal::IfNotEmpty;
      }
      throw IncorrectTypeInAssignment("The value for key \"" +
                                      std::string(key_) + "\" should be " +
                                      "\"Yes\", \"No\" or \"IfNotEmpty\".");
    }
  };

  /**
   * Flag to mark initialization with a YAML formatted string.
   */
  static const char InitializeFromYAMLString = 'S';

  /**
   * Flag to tune method(s) behavior such that it is descriptive from the
   * caller side. For example, see \ref extract_sub_configuration.
   */
  enum class GetEmpty { Yes, No };

  /**
   * Return type of Configuration::validate which conveys more information that
   * simply a two-state boolean variable.
   */
  enum class Is { Invalid, Deprecated, Valid };

  /**
   * Read config.yaml from the specified path.
   *
   * \param[in] path The directory where the SMASH config files are located.
   */
  explicit Configuration(const std::filesystem::path &path);

  /**
   * Read a YAML config file from the specified path.
   *
   * \param[in] path The directory where the SMASH config files are located.
   * \param[in] filename The filename (without path) of the YAML config file, in
   *                 case you don't want the default "config.yaml".
   */
  explicit Configuration(const std::filesystem::path &path,
                         const std::filesystem::path &filename);

  /**
   * Initialize configuration with a YAML formatted string.  This is
   * useful in 3-rd party application where we may not be able or
   * willing to read in external files.
   *
   * \param[in] yaml YAML formatted configuration data.
   * \param[in] sflag control flag InitializeFromYAMLString.
   */
  explicit Configuration(const char *yaml, const char sflag) {
    if (sflag == InitializeFromYAMLString) {
      merge_yaml(yaml);
    } else {
      throw std::runtime_error(
          "Unknown control flag in Configuration constructor"
          " with a YAML formatted string. Please, use"
          " Configuration::InitializeFromYAMLString.");
    }
  }

  /**
   * Trivial constructor for testing purposes
   *
   * \param[in] invalue Integer passed by value, currently it does nothing
   */
  explicit Configuration(int invalue) {
    if (invalue == 0) {
      invalue = 1;
    }
  }

#ifdef BUILD_TESTS
  /**
   * \mocking
   * Unit tests can use this constructor to get a Configuration object from a
   * built-in string.
   * This function is only available to tests and should never be used/needed in
   * actual SMASH code. The intention is to avoid creating a mock object for
   * Configuration to test other classes of SMASH.
   */
  explicit Configuration(const char *yaml) : root_node_(YAML::Load(yaml)) {}
#endif

  /// If you want to copy this you're doing it wrong
  Configuration(const Configuration &) = delete;
  /// If you want to copy this you're doing it wrong
  Configuration &operator=(const Configuration &) = delete;

  /// Moving is fine
  Configuration(Configuration &&) = default;
  /// Moving is fine
  Configuration &operator=(Configuration &&) = default;

  /**
   * Merge the configuration in \p yaml into the existing tree.
   *
   * The function parses the string in \p yaml into its internal tree
   * representation. Then it merges the nodes from the new tree into the
   * existing tree.
   * The merge resolves conflicts by taking the value from \p yaml.
   *
   * \param[in] yaml A string with YAML (or JSON) content that is to be merged.
   */
  void merge_yaml(const std::string &yaml);

  /// Lists all YAML::Nodes from the configuration setup.
  std::vector<std::string> list_upmost_nodes();

  /**
   * The default interface for SMASH to read configuration values.
   *
   * The function returns the value at the specified \p keys and removes it from
   * the Configuration object. Therefore, a subsequent call to take or has_value
   * with the same \p keys returns an undefined value / \c false.
   * By removing the value, the Configuration object keeps track which settings
   * were never read.
   *
   * \note If taking a key leaves the parent key without a value, then this is
   *       in turn removed and so on. From a performance point of view, it might
   *       be argued that this is not needed to be checkend and done at every
   *       \c take operation and it might be done once for all. However, on one
   *       hand it is a natural behaviour to expect and on the other hand this
   *       is hardly going to be an application bottle-neck.
   *
   * \param[in] keys You can pass an arbitrary number of keys inside curly
   * braces, following the nesting structure in the config file.
   * For example, given
   \verbatim
     Group:
         Key: Value
   \endverbatim
   * call \code string value = config.take({"Group", "Key"}); \endcode
   * to read the value. This will make the key \c "Group" also be removed from
   * the configuration, since it remains without any value.
   *
   * \return A proxy object that converts to the correct type automatically on
   *         assignment.
   */
  Value take(std::initializer_list<const char *> keys);

  /// \see take
  template <typename T>
  T take(std::initializer_list<const char *> keys, T default_value) {
    if (has_value(keys)) {
      return take(keys);
    }
    return default_value;
  }

  /**
   * Additional interface for SMASH to read configuration values without
   * removing them.
   *
   * The function returns the value at the specified \p keys but does not remove
   * it from the Configuration object. Semantically, this means the value was
   * not used.
   *
   * \param[in] keys You can pass an arbitrary number of keys inside curly
   * braces, following the nesting structure in the config file.
   *
   * \return A proxy object that converts to the correct type automatically on
   *         assignment.
   */
  Value read(std::initializer_list<const char *> keys) const;

  /// \see read
  template <typename T>
  T read(std::initializer_list<const char *> keys, T default_value) {
    if (has_value(keys)) {
      return read(keys);
    }
    return default_value;
  }

  /**
   * Remove all entries in the given section except for \p key.
   *
   * \param[in] key The key of the map entry to keep.
   * \param[in] section You can pass an arbitrary number of keys inside curly
   *                    braces, following the nesting structure in the config
   *                    file, in order to specify the section where to delete
   *                    entries. Omitting the \c section is equivalent to
   *                    specifying \c {} and the top-level section is
   *                    understood.
   */
  void remove_all_entries_in_section_but_one(
      const std::string &key, std::initializer_list<const char *> section = {});

  /**
   * Create a new configuration from a then-removed section of the present
   * object. This method is meant to be used to deal with sections only, i.e.
   * it will throw if used to extract a key value that is not a section (namely
   * a map in YAML language). Use \ref take for that purpose, instead.
   *
   * \param[in] keys You can pass an arbitrary number of keys inside curly
   *             braces, following the nesting structure in the config file.
   * \param[in] empty_if_not_existing
   *            Specify \c Configuration::GetEmpty::Yes if you want an empty
   *            Configuration in case the requested section does not exist.
   *
   * \throw std::runtime_error if the method is used
   *        - to access a scalar or sequence value;
   *        - to access a key that has no value or is an empty map;
   *        - to access a not existing key (unless explicitly allowed).
   *
   * \return A new \c Configuration containing the chosen section.
   */
  Configuration extract_sub_configuration(
      std::initializer_list<const char *> keys,
      Configuration::GetEmpty empty_if_not_existing =
          Configuration::GetEmpty::No);

  /**
   * Overwrite the value of the specified YAML node.
   *
   * \param[in] keys You can pass an arbitrary number of keys inside curly
   *                 braces, following the nesting structure in the config file.
   * \param[in] value An arbitrary value that yaml-cpp can convert into YAML
   *                  representation. Any builtin type, strings, maps, and
   *                  vectors can be used here.
   */
  template <typename T>
  void set_value(std::initializer_list<const char *> keys, T &&value) {
    auto node = find_node_at(root_node_, keys);
    node = std::forward<T>(value);
  }

  /**
   * Return whether there is a (maybe empty) value behind the requested \p keys.
   * \param[in] keys List of keys to be checked for
   */
  bool has_value_including_empty(
      std::initializer_list<const char *> keys) const;
  /**
   * Return whether there is a non-empty value behind the requested \p keys.
   * \param[in] keys List of keys to be checked for
   */
  bool has_value(std::initializer_list<const char *> keys) const;

  /**
   * @return \c true if the object is empty;
   * @return \c false if at least one key exists.
   */
  bool is_empty() const { return root_node_.size() == 0; }

  /**
   * Return a \c string of the current YAML tree.
   */
  std::string to_string() const;

  /**
   * Validate content of configuration in terms of YAML keys.
   *
   * A warning or error message is printed for deprecated or invalid keys,
   * respectively, together with information about SMASH versions, if possible.
   *
   * \note Here a full validation is done by default and all keys are checked,
   * although the validation might be shortened by returning \c false as soon as
   * an invalid key is found. However, a full validation is more user-friendly,
   * since as much information as possible about the input file is provided.
   *
   * \param[in] full_validation Whether all keys are checked or not.
   *
   * \return \c Is::Valid if the object contains valid keys only;
   * \return \c Is::Deprecated if the object is valid but has deprecated key(s);
   * \return \c Is::Invalid if the object contains at least one invalid key.
   */
  Is validate(bool full_validation = true) const;

 private:
  /** Create a subobject that has its root node at the given node.
   *
   * \note This constructor is not explicit because it can be called only from
   * inside Configuration and by making it explicit a return would require the
   * copy constructor.
   */
  Configuration(const YAML::Node &node)  // NOLINT(runtime/explicit) : see above
      : root_node_(YAML::Clone(node)) {}

  /**
   * Descend in the YAML tree from the given node using the provided keys.
   *
   * More precisely, it finds a node, copies its structure and replaces the
   * previous keys by the newly provided keys.
   *
   * \param[in] node YAML::Node to start the search from.
   * \param[in] keys Keys that will be used to descend the YAML tree.
   *
   * \return Node in the tree reached by using the provided keys.
   */
  YAML::Node find_node_at(YAML::Node node,
                          std::vector<const char *> keys) const;

  /// the general_config.yaml contents - fully parsed
  YAML::Node root_node_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_CONFIGURATION_H_
