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

#include <cassert>
#include <functional>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "stringfunctions.h"
#include "traits.h"

namespace smash {

/**
 * Descriptive alias for storing a SMASH version associated to keys metadata.
 * At the moment simply a \c std::string .
 */
using Version = std::string;

/**
 * Descriptive alias for storing keys metadata. At the moment this is only a
 * list of versions.
 */
using KeyMetadata = std::initializer_list<std::string_view>;

/**
 * Descriptive alias for storing key labels, i.e. the series of strings that
 * identify a key in the input file from the main section.
 * At the moment simply a \c std::vector<std::string> .
 */
using KeyLabels = std::vector<std::string>;

/**
 * Overload of the \c + operator to add a label to a \c KeyLabels object.
 *
 * @param lhs Constant lvalue reference to the key labels.
 * @param rhs Label to be added.
 * @return A new \c KeyLabels containing the desired labels.
 */
inline KeyLabels operator+(const KeyLabels& lhs, std::string_view rhs) {
  if (lhs.empty()) {
    return KeyLabels{std::string{rhs}};
  } else {
    KeyLabels result{lhs};
    result.push_back(std::string{rhs});
    return result;
  }
}

/**
 * Overload of the \c + operator to add a label to a \c KeyLabels object.
 *
 * @param lhs rvalue reference to the key labels.
 * @param rhs Label to be added.
 * @return A new \c KeyLabels containing the desired labels.
 */
inline KeyLabels operator+(KeyLabels&& lhs, std::string_view rhs) {
  if (lhs.empty()) {
    return KeyLabels{std::string{rhs}};
  } else {
    lhs.push_back(std::string{rhs});
    return std::move(lhs);  // See https://stackoverflow.com/a/14857144
  }
}

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
 * @tparam default_type Type of the key value. This \b must be a plain type, by
 *         that meaning have no cv-qualifier and not being any among the
 *         following types: array, pointer, function, or a mix of them.
 */
template <typename default_type>
class Key {
  static_assert(!std::is_const_v<default_type>);
  static_assert(!std::is_volatile_v<default_type>);
  static_assert(!std::is_reference_v<default_type>);
  static_assert(!std::is_pointer_v<default_type>);
  static_assert(!std::is_array_v<default_type>);
  static_assert(!std::is_function_v<default_type>);
  static_assert(!std::is_member_object_pointer_v<default_type>);
  static_assert(!std::is_member_function_pointer_v<default_type>);

 public:
  /**
   * \ingroup exception
   * Thrown when too few or too many versions are passed to the constructor.
   */
  struct WrongNumberOfVersions : public std::runtime_error {
    using std::runtime_error::runtime_error;
  };

  /**
   * @brief Descriptive alias for the key validator.
   */
  using validator_type = std::function<bool(const default_type&)>;

  /**
   * @brief Construct a new \c Key object without default value.
   *
   * @param[in] labels The label(s) identifying the key in the YAML input file.
   * @param[in] versions A list of one, two or three version numbers identifying
   * the versions in which the key has been introduced, deprecated and removed,
   * respectively.
   * @param[in] validator An optional functor that takes a default_type
   * parameters and returns a bool variable.
   *
   * @throw WrongNumberOfVersions If \c versions has the wrong size.
   */
  explicit Key(const KeyLabels& labels, const KeyMetadata& versions,
               validator_type validator = get_default_validator_())
      : Key{labels, Default<default_type>{}, versions, validator} {}

  /**
   * @brief Construct a new \c Key object with default value.
   *
   * @param[in] labels The label(s) identifying the key in the YAML input file.
   * @param[in] value The key default value.
   * @param[in] versions A list of one, two or three version numbers identifying
   * the versions in which the key has been introduced, deprecated and removed,
   * respectively.
   * @param[in] validator An optional functor that takes a default_type
   * parameters and returns a bool variable.
   *
   * @throw WrongNumberOfVersions If \c versions has the wrong size.
   * @throw std::invalid_argument If \c validator(value) returns \c false .
   */
  Key(const KeyLabels& labels, default_type value, const KeyMetadata& versions,
      validator_type validator = get_default_validator_())
      : Key{labels, Default<default_type>{value}, versions, validator} {}

  /**
   * @brief Construct a new \c Key object which is supposed to have a default
   * value, which however depends on other keys and will remain unset.
   *
   * @param[in] labels The label(s) identifying the key in the YAML input file.
   * @param[in] type_of_default The type of default value.
   * @param[in] versions A list of one, two or three version numbers identifying
   * the versions in which the key has been introduced, deprecated and removed,
   * respectively.
   * @param[in] validator An optional functor that takes a default_type
   * parameters and returns a bool variable.
   *
   * @throw WrongNumberOfVersions If \c versions has the wrong size.
   * @throw std::logic_error If \c type is not \c DefaultType::Dependent .
   */
  Key(const KeyLabels& labels, DefaultType type_of_default,
      const KeyMetadata& versions,
      validator_type validator = get_default_validator_())
      : Key{labels, Default<default_type>{type_of_default}, versions,
            validator} {}

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
   * @brief Get whether the given key value is valid.
   *
   * @param[in] value The value to be validated.
   *
   * @return \c true if the given value is valid,
   * @return \c false otherwise.
   */
  bool validate(const default_type& value) const noexcept {
    return validator_(value);
  }

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
   * Build and return a YAML-formatted string in the compact form (using braces
   * as single line).
   *
   * \param[in] value An \c std::optional value of the Key type. If a value is
   *            passed, this is added to the resulting string if its type is
   *            streamable using the \c << operator. If no value is passed and
   *            the key has a streamable default, this is used.
   *
   * @return \c std::string with labels formatted in a compact YAML format.
   */
  std::string as_yaml(
      std::optional<default_type> value = std::nullopt) const noexcept {
    std::stringstream value_as_string{};
    if constexpr (is_writable_to_stream_v<std::stringstream, default_type>) {
      if (value) {
        value_as_string << *value;
      } else if (default_.type_ == DefaultType::Value) {
        value_as_string << default_value();
      }
    }
    return as_yaml(value_as_string.str());
  }

  /**
   * Overload of the method taking a string as value. This can be useful for non
   * streamable types e.g. in tests.
   *
   * \note The passed \c value is not quoted and it is responsibility of the
   *       caller to properly quote it, if needed. This enables setting e.g.
   *       YAML maps as value.
   *
   * \see as_yaml
   */
  std::string as_yaml(std::string value) const noexcept {
    std::stringstream result{};
    result << "{" << smash::join(labels_, ": {") << ": " << value
           << smash::join(std::vector<std::string>(labels_.size(), "}"), "");
    return result.str();
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
    // Make nested class friend of enclosing one. This class is anyhow an
    // implementation detail and part of Key.
    friend class Key<T>;
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
  Key(const KeyLabels& labels, Default<default_type> value,
      const KeyMetadata& versions, validator_type validator)
      : default_{std::move(value)},
        labels_{labels.begin(), labels.end()},
        validator_{std::move(validator)} {
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
    /* Assert validator_ is set, which is usually the case unless in
     * particularly nasty scenarios to debug (e.g. calling this constructor from
     * a static/global object using a static/global validator and hence hitting
     * the undefined order of static initialisation). */
    assert(validator_);
    if (default_.value_ && validator_(*(default_.value_)) == false) {
      throw std::logic_error(
          "Key " + static_cast<std::string>(*this) +
          " has been declared with an invalid default value.");
    }
  }

  /**
   * Static method to get the default trivial validator.
   *
   * @return A functor that always returns \c true .
   *
   * @attention It might look unnecessary to have a static method and you might
   * think that the functor as a static member would be enough. However, this
   * would be in general wrong because this functor is used in the keys
   * constructors which are used by the \c InputKeys class, that is a collection
   * of static <tt>Keys</tt>. Hence, we would have static members using each
   * others and initialization of static members in C++ is undefined. We use
   * therefore the "construct on first use idiom", making the functor a static
   * object in a function scope. For more information, refer to <a
   * href="https://isocpp.org/wiki/faq/ctors#static-init-order-on-first-use-members">ISO
   * C++ FAQ</a>.
   */
  static validator_type get_default_validator_() {
    static auto always_true_functor = [](const default_type&) { return true; };
    return always_true_functor;
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
  /// The functor to validate key values
  validator_type validator_{get_default_validator_()};
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_KEY_H_
