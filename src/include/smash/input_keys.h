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

/*!\Userguide
 * \page validation Validation
 *
 */

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
 *            keys, you need to add it to the key_references_variant alias. In
 *            particular, you need to add a type to the \c std::variant which
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
 */
struct InputKeys {
  /*!\Userguide
   * \page validation Validation
   * \section General
   *
   * \anchor Modus \key Modus (string, required): \n
   * Selects a modus for the calculation, e.g.\ infinite matter
   * calculation, collision of two particles or collision of nuclei. The modus
   * will be configured in \ref input_modi_. Recognized values are:
   * \li \key Collider - For collisions of nuclei or compound objects.
   *     See \ref \ColliderModus
   * \li \key Sphere - For calculations of the expansion of a thermalized
   *      sphere. See \ref \SphereModus
   * \li \key Box - For infinite matter calculation in a rectangular box.
   *      See \ref \BoxModus
   * \li \key List - For given external particle list. See \ref \ListModus
   * \li \key ListBox - For given external particle list in the Box.
   */
  /**
   * See \ref Modus "user guide description" for more information.
   */
  inline static const Key<std::string> gen_modus{{"General", "Modus"}, {"1.0"}};
  /*!\Userguide
   * \page validation Validation
   * \anchor Delta_Time \key Delta_Time (double, optional, default: 1.0): \n
   * Fixed time step at which the collision-finding grid is recreated, and, if
   * potentials are on, momenta are updated according to the equations of
   * motion. The collision-finding grid finds all the collisions from time
   * t_{beginning_of_timestep} until time t_{beginning_of_timestep} +
   * Delta_Time, and puts them into a vector. The collisions are then sorted in
   * order of occurrence, and particles are propagated from collision to
   * collision. After each performed collision, additional collisions are found
   * for outgoing particles and merged into the sorted vector.
   */
  /**
   * See \ref Delta_Time "user guide description" for more information.
   */
  inline static const Key<double> gen_deltaTime{
      {"General", "Delta_Time"}, 1.0, {"1.0"}};

  /// Alias for the type to be used in the list of keys.
  using key_references_variant =
      std::variant<std::reference_wrapper<const Key<double>>,
                   std::reference_wrapper<const Key<std::string>>>;

  /// List of references to all existing SMASH keys.
  inline static const std::vector<key_references_variant> list = {
      std::cref(gen_modus), std::cref(gen_deltaTime)};
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_VALIDATION_H_