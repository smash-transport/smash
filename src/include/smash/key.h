/*
 *
 *    Copyright (c) 2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_KEY_H_
#define SRC_INCLUDE_SMASH_KEY_H_

#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "stringfunctions.h"

namespace smash {

/**
 * Descriptive alias for storing SMASH versions associated to keys metadata.
 * At the moment simply a \c std::string .
 */
using Version = std::string;

/**
 * Descriptive alias for storing key labels, i.e. the series of strings that
 * identify a key in the input file from the main section.
 * At the moment simply a \c std::vector<std::string> .
 */
using KeyLabels = std::vector<std::string>;

/**
 * @brief New type to explicit distinguish between mandatory and optional keys.
 */
enum class DefaultType {
  /// %Default "type" for mandatory keys
  Null,
  /// Normal default with a value associated to it
  Value,
  /// %Default value which depends on other keys
  Dependent
};

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
      : Key{labels, Default<default_type>{}, versions} {}

  /**
   * @brief Construct a new \c Key object with default value.
   *
   * @param[in] labels The label(s) identifying the key in the YAML input file.
   * @param[in] value The key default value.
   * @param[in] versions A list of one, two or three version numbers identifying
   * the versions in which the key has been introduced, deprecated and removed,
   * respectively.
   *
   * @throw WrongNumberOfVersions If \c versions has the wrong size.
   */
  Key(const std::initializer_list<std::string_view>& labels, default_type value,
      const std::initializer_list<std::string_view>& versions)
      : Key{labels, Default<default_type>{value}, versions} {}

  /**
   * @brief Construct a new \c Key object which is supposed to have a default
   * value, which however depends on other keys and will remain unset.
   *
   * @param[in] labels The label(s) identifying the key in the YAML input file.
   * @param[in] type_of_default The type of default value.
   * @param[in] versions A list of one, two or three version numbers identifying
   * the versions in which the key has been introduced, deprecated and removed,
   * respectively.
   *
   * @throw WrongNumberOfVersions If \c versions has the wrong size.
   * @throw std::logic_error If \c type is not \c DefaultType::Dependent .
   */
  Key(const std::initializer_list<std::string_view>& labels,
      DefaultType type_of_default,
      const std::initializer_list<std::string_view>& versions)
      : Key{labels, Default<default_type>{type_of_default}, versions} {}

  /**
   * @brief Let the clients of this class have access to the key type.
   */
  using type = default_type;

  /**
   * @brief Get the default value of the key.
   *
   * @return A \c default_type variable.
   *
   * @throw std::bad_optional_access If the key has no default value.
   */
  default_type default_value() const { return default_.value(); }

  /**
   * @brief Ask whether the default value depends on other other keys.
   *
   * @return \c true if this is the case,
   * @return \c false if the default value is known or the key is mandatory.
   */
  bool has_dependent_default() const noexcept {
    return default_.is_dependent();
  }

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
  bool has_same_labels(const KeyLabels& labels) const noexcept {
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

  /**
   * \brief Method to access the \c Key labels.
   *
   * \return A constant reference to the labels member for read-only access.
   */
  const KeyLabels& labels() const { return labels_; }

 private:
  /**
   * @brief Wrapper class around a type with the capability to both store the
   * type of default and its value, if any exists. This class has 3 valid
   * states:
   *
   * | State | `type_` | `value_` |
   * | :---: | :-----: | :------: |
   * | Required key | `DefaultType::Null`      | `std::nullopt`   |
   * | %Default value | `DefaultType::Value`     | ≠ `std::nullopt` |
   * | %Key with dependent default | `DefaultType::Dependent` | `std::nullopt` |
   *
   * There is a constructor for each of the cases above.
   *
   * @tparam T The default value type.
   *
   * \note This is an implementation detail of the \c Key class and it is meant
   * to be rigid in its usage. E.g., the constructor specifying a \c DefaultType
   * is meant to only accept \c DefaultType::Dependent because this is the only
   * way we want it to be used.
   */
  template <typename T>
  class Default {
   public:
    /**
     * @brief Construct a new \c Default object which denotes a mandatory value
     * without a default. This is meant to be used for required keys.
     */
    Default() : type_{DefaultType::Null} {}
    /**
     * @brief Construct a new \c Default object storing its default value.
     *
     * @param in The default value to be stored
     */
    explicit Default(T in) : value_{std::move(in)} {}
    /**
     * @brief Construct a new \c Default object which has a value dependent on
     * external information.
     *
     * @param type The type of default (it should be \c DefaultType::Dependent
     * ).
     *
     * @throw std::logic_error if called with a type different from \c
     * DefaultType::Dependent .
     */
    explicit Default(DefaultType type) : type_{type} {
      if (type != DefaultType::Dependent) {
        throw std::logic_error("Default constructor used with invalid type!");
      }
    }

    /**
     * @brief Retrieve the default value stored in the object
     *
     * @return The default value stored
     *
     * @throw std::bad_optional_access If the object stores no default value.
     */
    T value() const { return value_.value(); }

    /**
     * @brief Ask whether the default value depends on other external
     * information.
     *
     * @return \c true if this is the case,
     * @return \c false if the default value is known or none exists.
     */
    bool is_dependent() const noexcept {
      return type_ == DefaultType::Dependent;
    }

   private:
    /// The type of default value
    DefaultType type_ = DefaultType::Value;
    /// The default value, if any
    std::optional<T> value_ = std::nullopt;
  };

  /**
   * @brief Private constructor of the Key object.
   *
   * This is meant to do the real construction, while the other public
   * constructors just delegate to this one. This is possible because this
   * constructor takes a \c Default argument and the other construct one to
   * delegate construction.
   *
   * @see public constructor documentation for the parameters description.
   */
  Key(const std::initializer_list<std::string_view>& labels,
      Default<default_type> value,
      const std::initializer_list<std::string_view>& versions)
      : default_{std::move(value)}, labels_{labels.begin(), labels.end()} {
    /*
     * The following switch statement is a compact way to initialize the
     * three version member variables without repetition and lots of logic
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

  /// SMASH version in which the key has been introduced
  Version introduced_in_{};
  /// SMASH version in which the key has been deprecated, if any
  std::optional<Version> deprecated_in_{};
  /// SMASH version in which the key has been removed, if any
  std::optional<Version> removed_in_{};
  /// Key default value
  Default<default_type> default_{};
  /// The label(s) identifying the key in the YAML input file
  KeyLabels labels_{};
};

}

#endif