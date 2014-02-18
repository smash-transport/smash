/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_CONFIGURATION_H_
#define SRC_INCLUDE_CONFIGURATION_H_

#include <stdexcept>
#include <string>

#include <yaml-cpp/yaml.h>

#ifndef DOXYGEN
namespace boost {
namespace filesystem {
class path;
}  // namespace filesystem
}  // namespace boost
#endif

/**
 * Interface to the SMASH configuration files.
 */
class Configuration {
 public:
  /// Thrown when the types in the config file and C++ don't match.
  struct IncorrectTypeInAssignment : public std::runtime_error {
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
    /// a YAML leaf node
    const YAML::Node node_;

   public:
    Value(const YAML::Node &n) : node_(n) {
      assert(n.IsScalar() || n.IsSequence() || n.IsMap());
    }

    /// if you want to copy this you're doing it wrong
    Value(const Value &) = delete;
    /// if you want to copy this you're doing it wrong
    Value &operator=(const Value &) = delete;

    /**
     * This function determines the type it is assigned to and calls
     * YAML::Node::as<T>() with this type.
     *
     * This makes reading values more convenient than calling as<type>()
     * explicitly.
     */
    template <typename T>
    operator T() {
      try {
        return node_.as<T>();
      }
      catch (YAML::TypedBadConversion<T> &e) {
        throw IncorrectTypeInAssignment(
            "The value cannot be converted to the requested type.");
      }
    }
  };

  /**
   * Reads config_general.yaml from the specified \p path.
   *
   * \param path The directory where the SMASH config files are located.
   */
  explicit Configuration(const boost::filesystem::path &path);

  /**
   * Reads a YAML config file from the specified \p path.
   *
   * \param path The directory where the SMASH config files are located.
   * \param filename The filename (without path) of the YAML config file, in
   *                 case you don't want the default "config_general.yaml".
   */
  explicit Configuration(const boost::filesystem::path &path,
                         const std::string &filename);

  Configuration(const Configuration &) = delete;
  Configuration &operator=(const Configuration &) = delete;

  Configuration(Configuration &&) = default;
  Configuration &operator=(Configuration &&) = default;

  void merge_yaml(const std::string &yaml);

  /**
   * The default interface for SMASH to read configuration values.
   *
   * The function returns the value at the specified \p keys and removes it from
   * the Configuration object. Therefore, a subsequent call to take or has_value
   * with the same \p keys returns an undefined value / \c false.
   * By removing the value, the Configuration object keeps track which settings
   * were never read.
   *
   * \param keys You can pass an arbitrary number of keys inside curly braces,
   *             following the nesting structure in the config file. Example:
                 \verbatim
     Group:
         Key: Value
                 \endverbatim
   *             Call \code string value = config.take({"Group", "Key"});
   *             \endcode to read the value.
   *
   * \return A proxy object that converts to the correct type automatically on
   *         assignment.
   */
  Value take(std::initializer_list<const char *> keys);

  /**
   * Additional interface for SMASH to read configuration values without
   * removing them.
   *
   * The function returns the value at the specified \p keys but does not remove it from
   * the Configuration object. Semantically, this means the value was not used.
   *
   * \param keys You can pass an arbitrary number of keys inside curly braces,
   *             following the nesting structure in the config file.
   *
   * \return A proxy object that converts to the correct type automatically on
   *         assignment.
   */
  Value read(std::initializer_list<const char *> keys) const;

  /**
   * Access to the YAML::Node behind the requested \p keys.
   *
   * If you want to read a value use the \ref read function above. Use the
   * subscript operator if you want to assign a new value. The YAML::Node class
   * will automatically convert the data you assign to a string representation
   * suitable for the YAML file.
   *
   * \return An opaque object that can be assigned to.
   *
   * \see take
   * \see read
   */
  template <typename T>
  Configuration operator[](T &&key) {
    return root_node_[std::forward<T>(key)];
  }

  template <typename T>
  Configuration &operator=(T &&value) {
    root_node_ = std::forward<T>(value);
    return *this;
  }

  /**
   * Returns whether there is a value behind the requested \p keys.
   */
  bool has_value(std::initializer_list<const char *> keys) const;

  /**
   * Returns a string listing the key/value pairs that have not been taken yet.
   */
  std::string unused_values_report() const;

 private:
  /// Creates a subobject that has its root node at the given node.
  Configuration(const YAML::Node &node) : root_node_(node) {
  }

  /// the general_config.yaml contents - fully parsed
  YAML::Node root_node_;
};
#endif  // SRC_INCLUDE_CONFIGURATION_H_
