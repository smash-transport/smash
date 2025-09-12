/*
 *
 *    Copyright (c) 2014-2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_CONFIGURATION_H_
#define SRC_INCLUDE_SMASH_CONFIGURATION_H_

#include <algorithm>
#include <array>
#include <exception>
#include <filesystem>
#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "yaml-cpp/yaml.h"

#include "cxx17compat.h"
#include "forwarddeclarations.h"
#include "key.h"
#include "stringify.h"

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
  static Node encode(const T &x) {
    if constexpr (std::is_convertible_v<T, std::string>) {
      return Node{static_cast<std::string>(x)};
    } else {
      static_assert(smash::has_to_string_v<T>,
                    "Encoding type T to YAML::Node requires an overload of "
                    "smash::to_string(T) to convert T to an std::string.");
      return Node{smash::to_string(x)};
    }
  }

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
 * \page doxypage_input_particles
 *
 * <h3>How the particles file is used</h3>
 * The particles available to SMASH are defined in the *input/particles.txt*
 * file. The content of this file is internally copied by CMake to the
 * ***build*** directory when running the `cmake` command **for the first time**
 * to set up SMASH. If you want to modify the particles file, you are encouraged
 * to copy the provided one to a wished location, which has then to be passed to
 * SMASH via the `-p` option. For example, assuming to have a
 * *custom_particles.txt* file in the ***build*** folder, the SMASH executable
 * can be run from there and instructed to use the own particles file via
 * ```console
 * ./smash -p custom_particles.txt
 * ```
 *
 * <h3>The particle file format</h3>
 * %Particles are specified as a table with particles properties in different
 * columns, which may be separated by an arbitrary number of spaces:
 * ```
 * <name> <mass in GeV> <width in GeV> <parity> <PDG codes>
 * ```
 * The name has to be a unique UTF-8 string. Conventionally, unicode names are
 * used in SMASH to make the file more readable and generate prettier output. It
 * is possible to only specify the isospin multiplet and SMASH will fill in the
 * properties of the components of the multiplet assuming isospin symmetry. The
 * names generated this way will have the charges appended to the multiplet name
 * using the unicode characters `⁻`, `⁰` and `⁺`. This is appropriate for almost
 * all particles. Anti particles do not have to be specified explicitly.
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
 * specify the π multiplet using the following line in *particles.txt*, where
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
 * Comments can be added to the particles file using the `#` character.
 * Everything after `#` until the end of the line is ignored.
 *
 * <hr>
 * \attention
 * -# If you specify an incorrect value, SMASH will print an error similar to
 *    the following:
 *    ```
 *    Failed to convert the input string to the expected data types.
 *    ```
 * -# SMASH validates (up to some small numeric precision) the mass of some
 *    particles (e.g. N, π, K, ω, Δ, D). Therefore, totally nonphysical mass
 *    values cannot be used and SMASH will abort with a message error like e.g.
 *    the following
 *    ```
 *    Nucleon mass in input file different from 0.938000
 *    ```
 *    if a mass discrepancy for one of the validated particle is detected. This
 *    constraint is due to internal assumptions in SMASH and is needed to
 *    guarantee the correctness of the results. If you really need to use SMASH
 *    with nonphysical mass values, feel free to contact us or open an issue.
 * -# Related to the previous point, it is important to mention that all hadrons
 *    belonging to the same isospin multiplet must have the same mass and this
 *    is enforced by SMASH, which will fail otherwise. Feel free to get in touch
 *    with us, if this restriction represents a problem for you.
 * -# Some reactions in SMASH are parametrized and require specific particles in
 *    the final state. When such a reaction happens and the required particle is
 *    not defined, SMASH will crash.
 * -# When running a box simulation in which detailed balance is expected to be
 *    conserved, the particles file will need to be modified. See \ref
 *    modi_box_usage_remark "this remark about the box modus" for further
 *    information.
 */

/*!\Userguide
 * \page doxypage_input_decaymodes
 *
 * All possible decays and resonance formations in SMASH are provided by the
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
 * The names have to be the ones defined in *particles.txt* (see \ref
 * doxypage_input_particles). If multiplet names are used, the other branching
 * ratios are generated by SMASH assuming isospin symmetry. Note that currently
 * decay channels can only be specified for whole multiplets; individual
 * particles can however still be used in a decay channel as specific daughters.
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
 * \note
 * Note that SMASH has an internal width cutoff (currently 10 keV), below which
 * particles cannot decay, even if decays are specified in `decaymodes.txt`.
 * This only changes when `Ignore_Minimum_Decay_Width_For_Decays_At_The_End` is
 * enabled, then all decays are considered in the final decays, no matter their
 * decay width.
 *
 * \note
 * Note further, that the decay modes file will need to be modified when running
 * a box simulation in which detailed balance is expected to be conserved. See
 * \ref doxypage_input_conf_modi_box for further information.
 */

/**
 * Interface to the SMASH configuration files.
 *
 * The configuration is created from a %YAML file and then stores a nested map
 * of maps (normally a tree, but %YAML allows it to be cyclic - even though we
 * don't want that feature). Since the resource owned by the object is a
 * \c YAML::Node that handle memory in a similar way as a pointer does, it is
 * forbidden (nor should it be needed) to copy instances of this class, while
 * moving is fine (see special members documentation for more information).
 *
 * The typical usage of a Configuration is to create it, consume (i.e.
 * <tt>take</tt>) all its values and let it being destructed. Since this is the
 * contact point with SMASH input file, the class is meant to be strict in its
 * usage, so that it is possible to help the inexpert user, who might being
 * using a wrong input file and/or e.g. specify an unused key hoping in an
 * effect that indeed does not occur. Therefore, it is imposed that <b>all keys
 * must be parsed before an instance gets destroyed</b>. If this is not the
 * case, an exception will be thrown.
 *
 * For the typical usage in SMASH one needs to read the value once. In that
 * case, use the Configuration::take function, for example:
 * \code
 * Key<double> key{{"General", "SIGMA"}, {"1.0"}};
 * double value = config.take(key);
 * \endcode
 * Note that only Key objects can be taken. This is done on purpose, because all
 * SMASH allowed keys are gathered as static members of the InputKeys class.
 * This allows calls like \c config.take(InputKeys::key_name) (where of course
 * an existing member should be used).
 *
 * The opposite operation of \c take is the Configuration::set_value method,
 * which has a similar syntax, but needs the new value to be assigned, e.g.
 * \code
 * config.set_value(key, 3.1415);
 * \endcode
 *
 * If you need to delegate parsing of a section to some object, you can use the
 * Configuration::extract_complete_sub_configuration method, which is taking a
 * full section and returning a new, distinct Configuration instance.
 *
 * Last but not least, the Configuration::validate method is used by SMASH to
 * check that all given keys are allowed in the present version of the codebase.
 * This is achieved by querying the "database" InputKeys class.
 *
 * \attention As the Configuration is implemented, it does not make sense in
 * practice to have constant instances, because their keys could not be taken
 * and their destruction would lead to an exception being thrown. However, it
 * still makes perfectly sense to have constant methods (think e.g. of a
 * <tt>const %Configuration&</tt> being passed to a function).
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
   * \ingroup exception
   * Thrown if a Key is taken twice.
   */
  struct TakeSameKeyTwice : public std::logic_error {
    using std::logic_error::logic_error;
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

#ifdef BUILD_TESTS
  /**
   * \mocking
   * Unit tests can use this constructor to get a Configuration object from a
   * built-in string.
   * This function is only available to tests and should never be used/needed in
   * actual SMASH code. The intention is to avoid creating a mock object for
   * Configuration to test other classes of SMASH.
   */
  explicit Configuration(const char *yaml) : root_node_(YAML::Load(yaml)) {
    if (root_node_.IsNull())
      root_node_ = YAML::Node{YAML::NodeType::Map};
  }
#endif

  /**
   * Prevent Configuration objects from being copied.
   *
   * Underneath, the resource is a \c YAML::Node and since this handles memory
   * in a similar way as a pointer would do, copying an object would make
   * several instances point to the same memory and it would make it difficult
   * to use this object correctly. Therefore, copies are not allowed.
   */
  Configuration(const Configuration &) = delete;
  /**
   * Prevent Configuration objects from being copy-assigned.
   *
   * See copy constructor Configuration(const Configuration &) for more
   * information.
   */
  Configuration &operator=(const Configuration &) = delete;

  /**
   * Provide class with move constructor.
   *
   * In contrast to copying, moving is fine, since this keeps the owner
   * of the resource unique.
   *
   * \note Since the class has the peculiar behavior that all keys must be
   * parsed before it gets destroyed (otherwise an exception is thrown),
   * it is important to manually implement the move operations, in order to
   * ensure that objects that are moved from result cleared and their
   * destruction is not leading to any throw. This is not guaranteed if the
   * special members are defaulted to the compiler generated versions.
   */
  Configuration(Configuration &&);

  /**
   * Provide class with move assignment operator.
   *
   * See move constructor Configuration(Configuration &&) for more information.
   */
  Configuration &operator=(Configuration &&);

  /**
   * Destroy the object, optionally throwing if not all keys were taken.
   *
   * This is a way to enforce that the object has to be consumed (i.e.
   * completely parsed) during its lifetime. Since this object might be
   * destructed during stack unwinding, the destructor has to throw only
   * if it is safe to do so and the uncaught_exceptions_ member is used
   * to properly implement this behavior.
   */
  ~Configuration() noexcept(false);

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
   * The function returns the value at the specified \p Key and removes its
   * labels from the Configuration object. Therefore, a subsequent call to the
   * \c take or \c has_value methods with the same \p Key throws or returns
   * \c false respectively. By removing the value, the Configuration object
   * keeps track which keys were never taken.
   *
   * \attention If a not existent Key is taken, its default value is returned,
   *            if any exists. If a not existing required key is taken, an error
   *            will be given.
   *
   * \note If taking a key leaves the parent key without a value, then this is
   *       in turn removed and so on. From a performance point of view, it might
   *       be argued that this is not needed to be checked and done at every
   *       \c take operation and it might be done once for all. However, on one
   *       hand it is a natural behaviour to expect and on the other hand this
   *       is hardly going to be an application bottle-neck.
   *
   * \warning Since \c take returns the default value when the key is not
   *          present in the configuration, it is important to make it throw if
   *          an existing key is attempted to be taken twice. Otherwise it would
   *          happen that taking any existing key would return the user-defined
   *          value the first time and taking it again would return the key
   *          default value. This is a misleading behaviour we want to avoid.
   *          However, a key can be taken several times if it exists at take
   *          time. For example, taking a key, setting its value and taking it
   *          again is a valid behaviour.
   *
   * \param[in] key The input key that should be taken. This is usually one of
   * the \c InputKeys static members, i.e. one of the allowed keys. Of course,
   * any \c Key would work. For example, given
   \verbatim
   Group:
     Key: 42
   \endverbatim
   * then
   \verbatim
   Key key<int>{{"Group", "Key"}, {"1.0"}};
   string value = config.take(key);
   \endverbatim
   * will take the value. This will make the key \c "Group" also be removed from
   * the configuration, since it remains without any value.
   *
   * \return The value of the taken key if present, its default value otherwise.
   *
   * \throw TakeSameKeyTwice if a key was already previously taken.
   * \throw std::invalid_argument if a key without a default is taken but it is
   *        absent in the configuration.
   */
  template <typename T>
  T take(const Key<T> &key) {
    if (has_value(key)) {
      // The following return statement converts a Value into T
      return take({key.labels().begin(), key.labels().end()});
    } else if (has_section(key.labels())) {
      // In this case, if the Key type is a map, we take it, otherwise fails
      if constexpr (isMap<typename Key<T>::type>::value) {
        return take({key.labels().begin(), key.labels().end()});
      } else {
        throw std::logic_error(
            "Key " + std::string{key} +  // NOLINT(whitespace/braces)
            " was taken, but its value is not a map, although there is a "
            "section in the configuration with its labels.");
      }
    } else if (did_key_exist_and_was_it_already_taken(key.labels())) {
      throw TakeSameKeyTwice("Attempt to take key " +
                             std::string{key} +  // NOLINT(whitespace/braces)
                             " twice.");
    } else {
      try {
        return key.default_value();
      } catch (std::bad_optional_access &) {
        throw std::invalid_argument(
            "Key " + std::string{key} +  // NOLINT(whitespace/braces)
            " without default value taken, but missing in configuration.");
      }
    }
  }

  /**
   * Alternative method to take a key value, specifying the default value.
   * \see take
   *
   * @tparam T The type of the key to be taken
   * @param key The key to be taken
   * @param default_value The default value to be returned if the key is not
   *                      present in the configuration
   * @return The value of the key
   *
   * \throw std::logic_error If the key has not a default value declared as
   * dependent on external entities.
   */
  template <typename T>
  T take(const Key<T> &key, T default_value) {
    if (!key.has_dependent_default()) {
      throw std::logic_error(
          "An input Key without dependent default cannot be taken specifying a "
          "default value! Either define the key as having a dependent default "
          "or take it without a default value (which is a Key property).");
    }
    if (has_value(key)) {
      return take(key);
    }
    return default_value;
  }

  /**
   * Additional interface for SMASH to read configuration values without
   * removing them.
   *
   * The function returns the value of the specified \c Key but does not remove
   * it from the Configuration object. Semantically, this means the value was
   * not used.
   *
   * Also this method returns the default value of the key (or an error if none
   * is available), if the key is not present in the configuration. \see take
   *
   * \note Since reading a key does not remove it from the configuration, it is
   *       not necessary to store the key read to avoid taking it multiple
   *       times. Actually, doing so is safe and will return the same value.
   *
   * \param[in] key The input key that should be taken.
   *
   * \return The value of the taken key if present, its default value otherwise.
   *
   * \throw std::invalid_argument if a key without a default is taken but it is
   *        absent in the configuration.
   */
  template <typename T>
  T read(const Key<T> &key) const {
    if (has_value(key)) {
      // The following return statement converts a Value into T
      return read({key.labels().begin(), key.labels().end()});
    } else if (has_section(key.labels())) {
      // In this case, if the Key type is a map, we take it, otherwise fails
      if constexpr (isMap<typename Key<T>::type>::value) {
        return read({key.labels().begin(), key.labels().end()});
      } else {
        throw std::logic_error(
            "Key " + std::string{key} +  // NOLINT(whitespace/braces)
            " was read, but its value is not a map, although there is a "
            "section in the configuration with its labels.");
      }
    } else {
      try {
        return key.default_value();
      } catch (std::bad_optional_access &) {
        throw std::invalid_argument(
            "Key " + std::string{key} +  // NOLINT(whitespace/braces)
            " without default value read, but missing in configuration.");
      }
    }
  }

  /**
   * Alternative method to read a key value, specifying the default value.
   * \see read
   *
   * @tparam T The type of the key to be read
   * @param key The key to be read
   * @param default_value The default value to be returned if the key is not
   *                      present in the configuration
   * @return The value of the key
   *
   * \throw std::logic_error If the key has not a default value declared as
   * dependent on external entities.
   */
  template <typename T>
  T read(const Key<T> &key, T default_value) {
    if (!key.has_dependent_default()) {
      throw std::logic_error(
          "An input Key without dependent default cannot be read specifying a "
          "default value! Either define the key as having a dependent default "
          "or read it without a default value (which is a Key property).");
    }
    if (has_value(key)) {
      return read(key);
    }
    return default_value;
  }

  /**
   * Overwrite the value of the YAML node corresponding to the specified key.
   *
   * \param[in] key The input key that should be changed.
   * \param[in] value An arbitrary value that yaml-cpp can convert into YAML
   *                  representation. Any builtin type, strings, maps, and
   *                  vectors can be used here. Of course, this has to match the
   *                  Key type passed as first argument.
   *
   * \tparam T The type of the value to be assigned to the Key.
   * \tparam U The type of the key value. This is by default \c T but
   *           it has been allowed to be different from it, as long as it is
   *           convertible to T. This enables e.g. to set a key with a string
   *           value using a <tt>const char*</tt> second argument.
   *
   * \attention This method creates a new entry in the configuration if the
   *            passed key is not yet existing in it.
   *
   * \note Removing qualifiers and the reference in the default value of the
   *       second template argument is needed because the Key type is a plain
   *       type and \c T might be deduced to a constant and/or reference type.
   */
  template <typename T, typename U = remove_cvref_t<T>,
            typename std::enable_if_t<std::is_convertible_v<T, U>, bool> = true>
  void set_value(Key<U> key, T &&value) {
    auto node = find_node_creating_it_if_not_existing(
        {key.labels().begin(), key.labels().end()});
    node = std::forward<T>(value);
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
  void remove_all_entries_in_section_but_one(const std::string &key,
                                             KeyLabels section = {});

  /**
   * Create a new configuration from a then-removed section of the present
   * object. This method is meant to be used to deal with sections only, i.e.
   * it will throw if used to extract a key value that is not a section (namely
   * a map in YAML language). Use \ref take for that purpose, instead.
   *
   * \param[in] section You can pass an arbitrary number of keys inside curly
   *            braces, following the nesting structure in the config file.
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
      KeyLabels section, Configuration::GetEmpty empty_if_not_existing =
                             Configuration::GetEmpty::No);

  /**
   * Alternative method to extract a sub-configuration, which retains the labels
   * from the top-level in the returned object instead of dropping them.
   *
   * \see extract_sub_configuration
   */
  Configuration extract_complete_sub_configuration(
      KeyLabels section, Configuration::GetEmpty empty_if_not_existing =
                             Configuration::GetEmpty::No);

  /**
   * Enclose the configuration into the given section.
   *
   * @param section The section in which to enclose the configuration.
   */
  void enclose_into_section(KeyLabels section);

  /**
   * Return whether the configuration has a (possibly empty) non-map key.
   * Although %YAML keys can have maps as value, we rather refer to those as
   * sections and we do not consider them as key in the SMASH database sense. A
   * key has then either a scalar or sequence value.
   *
   * \param[in] key The key to be checked for.
   */
  template <typename T>
  bool has_key(const Key<T> &key) const {
    const auto found_node =
        find_existing_node({key.labels().begin(), key.labels().end()});
    return found_node.has_value() && !(found_node.value().IsMap());
  }

  /**
   * Return whether there is a \b non-empty value behind the requested \p key
   * (which is supposed not to refer to a section). If there is a section with
   * the same labels as the provided key has, this function returns \c false .
   *
   * \param[in] key The key to be checked for.
   */
  template <typename T>
  bool has_value(const Key<T> &key) const {
    const auto found_node =
        find_existing_node({key.labels().begin(), key.labels().end()});
    return found_node.has_value() && !(found_node.value().IsNull()) &&
           !(found_node.value().IsMap());
  }

  /**
   * Return whether there is a (possibly empty) section with the given labels.
   *
   * \param[in] labels The labels of the section to be checked for.
   */
  bool has_section(const KeyLabels &labels) const {
    const auto found_node = find_existing_node({labels.begin(), labels.end()});
    return found_node.has_value() && found_node.value().IsMap();
  }

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
   * Erase the Configuration content.
   *
   * This function is useful e.g. in tests to clean up not taken keys
   * that would trigger an exception being thrown at by the destructor.
   */
  void clear() { root_node_.reset(); }

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
  /**
   * Proxy object to be used when taking or reading keys in the configuration.
   * This type automatically converts to the target type e.g. on assignment. An
   * object of this type is returned by the private \c take and \c read methods
   * and it is constructed from the \c YAML::Node to be taken or read.
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
     * Set total cross section strategy from configuration values.
     *
     * \return TotalCrossSectionStrategy.
     * \throw IncorrectTypeInAssignment in case a strategy that is
     * not available is provided as a configuration value.
     */
    operator TotalCrossSectionStrategy() const {
      const std::string s = operator std::string();
      if (s == "BottomUp") {
        return TotalCrossSectionStrategy::BottomUp;
      }
      if (s == "TopDown") {
        return TotalCrossSectionStrategy::TopDown;
      }
      if (s == "TopDownMeasured") {
        return TotalCrossSectionStrategy::TopDownMeasured;
      }
      throw IncorrectTypeInAssignment(
          "The value for key \"" + std::string(key_) + "\" should be " +
          "\"BottomUp\", \"TopDown\" " + "or \"TopDownMeasured\".");
    }

    /**
     * Set how pseudo-resonances are used from configuration values.
     *
     * \return PseudoResonance.
     * \throw IncorrectTypeInAssignment in case a key that is
     * not available is provided as a configuration value.
     */
    operator PseudoResonance() const {
      const std::string s = operator std::string();
      if (s == "None") {
        return PseudoResonance::None;
      }
      if (s == "Largest") {
        return PseudoResonance::Largest;
      }
      if (s == "Closest") {
        return PseudoResonance::Closest;
      }
      if (s == "LargestFromUnstable") {
        return PseudoResonance::LargestFromUnstable;
      }
      if (s == "ClosestFromUnstable") {
        return PseudoResonance::ClosestFromUnstable;
      }
      throw IncorrectTypeInAssignment(
          "The value for key \"" + std::string(key_) + "\" should be " +
          "\"None\",  \"Largest\",  \"Closest\", \"LargestFromUnstable\", or "
          "\"ClosestFromUnstable\".");
    }

    /**
     * Set condition of fluidization for hydrodynamic initial
     * conditions.
     *
     * \return FluidizationType.
     * \throw IncorrectTypeInAssignment in case a key that is
     * not available is provided as a configuration value.
     */
    operator FluidizationType() const {
      const std::string s = operator std::string();
      if (s == "Constant_Tau") {
        return FluidizationType::ConstantTau;
      } else if (s == "Dynamic") {
        return FluidizationType::Dynamic;
      }
      throw IncorrectTypeInAssignment("The value for key \"" +
                                      std::string(key_) + "\" should be " +
                                      "\"Constant_Tau\" or \"Dynamic\".");
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

    /**
     * Set FluidizableProcessesBitSet from configuration values.
     *
     * \return FluidizableProcessesBitSet with all included reaction types.
     * \throw IncorrectTypeInAssignment in case a reaction type that is not
     * available is provided as a configuration value.
     */
    operator FluidizableProcessesBitSet() const {
      const std::vector<std::string> v = operator std::vector<std::string>();
      FluidizableProcessesBitSet s;
      for (const auto &x : v) {
        if (x == "All") {
          s.set();
          break;
        } else if (x == "Elastic") {
          s.set(IncludedFluidizableProcesses::From_Elastic);
        } else if (x == "Decay") {
          s.set(IncludedFluidizableProcesses::From_Decay);
        } else if (x == "Inelastic") {
          s.set(IncludedFluidizableProcesses::From_Inelastic);
        } else if (x == "SoftString") {
          s.set(IncludedFluidizableProcesses::From_SoftString);
        } else if (x == "HardString") {
          s.set(IncludedFluidizableProcesses::From_HardString);
        } else {
          throw IncorrectTypeInAssignment(
              "The value for key \"" + std::string(key_) +
              "\" should be \"All\", \"Elastic\", \"Decay\", "
              "\"Inelastic\", \"SoftString\", \"HardString\", "
              "or any combination of these.");
        }
      }
      return s;
    }
  };

  /**
   * Create a sub-object that has its root node at the given node.
   *
   * \note This constructor is not explicit because it can be called only from
   * inside Configuration and by making it explicit a return would require the
   * copy constructor.
   */
  Configuration(const YAML::Node &node)  // NOLINT(runtime/explicit) : see above
      : root_node_(YAML::Clone(node)) {}

  /**
   * Descend in and if needed modify the YAML tree from the given node using the
   * provided keys.
   *
   * After this call nodes corresponding to the passed keys are guaranteed to
   * exist in the tree.
   *
   * \param[in] keys Keys that will be possibly added to the YAML tree.
   *
   * \return Node in the tree reached by using the provided keys.
   */
  YAML::Node find_node_creating_it_if_not_existing(
      std::vector<std::string_view> keys) const;

  /**
   * Descend in the YAML tree from the given node using the provided keys.
   *
   * This function \b must not use the YAML::Node subscript operator, which is
   * at the very bottom level creating an undefined node in the %YAML tree,
   * hence "wasting" some memory. Note that the fact that this method is marked
   * as const does not forbid to use the access operator on <tt>root_node_</tt>,
   * because of how the %YAML library works. We want the tree to be completely
   * untouched by this method.
   *
   * \param[in] keys Keys that will be used to descend the YAML tree.
   *
   * \return \c std::optional<YAML::Node> containing the node in the tree
   *         reached by using the provided keys, if it exists;
   * \return \c std::nullopt otherwise.
   *
   * \note It has been decided to return an optional value rather than throwing
   *       an exception because this method is going to be used in other methods
   *       like \c has_value and putting there a try-catch block would probably
   *       cause a performance cost that can be avoided (exceptions on the
   *       exceptional path are expensive).
   */
  std::optional<YAML::Node> find_existing_node(
      std::vector<std::string_view> keys) const;

  /**
   * Utility type trait (general case) for the take and read public methods.
   */
  template <class T>
  struct isMap : std::false_type {};
  /**
   * Utility type trait (special case) for the take and read public methods.
   */
  template <class Key, class Value>
  struct isMap<std::map<Key, Value>> : std::true_type {};

  /**
   * This is the implementation detail to take a key. Having a non-templated
   * method for it it allows for defining the method in the source file, which
   * isn't possible for a template.
   *
   * \param[in] labels the labels of the key to be taken
   *
   * \return The \c Value of the key
   *
   * \throw std::runtime_error if the key does not exist
   */
  Value take(std::vector<std::string_view> labels);

  /**
   * This is the implementation detail to read a key. Having a non-templated
   * method for it it allows for defining the method in the source file, which
   * isn't possible for a template.
   *
   * \param[in] labels the labels of the key to be read
   *
   * \return The \c Value of the key
   *
   * \throw std::runtime_error if the key does not exist
   */
  Value read(std::vector<std::string_view> labels) const;

  /**
   * Find out whether a key has been already taken.
   *
   * \param labels The labels of the key to be checked
   * \return \c true if the key was already taken,
   * \return \c false otherwise.
   */
  bool did_key_exist_and_was_it_already_taken(const KeyLabels &labels) const {
    return std::find(existing_keys_already_taken_.begin(),
                     existing_keys_already_taken_.end(),
                     labels) != existing_keys_already_taken_.end();
  }

  /// The general_config.yaml contents - fully parsed
  YAML::Node root_node_{YAML::NodeType::Map};

  /// Counter to be able to optionally throw in destructor
  int uncaught_exceptions_{std::uncaught_exceptions()};

  /// List of taken keys to throw on taking same key twice
  std::vector<KeyLabels> existing_keys_already_taken_{};
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_CONFIGURATION_H_
