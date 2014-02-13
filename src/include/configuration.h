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

namespace boost {
namespace filesystem {
class path;
}  // namespace filesystem
}  // namespace boost

/**
 * Interface to the SMASH configuration files.
 */
class Configuration {
  /// the general_config.yaml contents - fully parsed
  YAML::Node config_general_;  // this needs to be on top for decltype in
                               // operator[] to compile.

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
    const YAML::Node &node_;

   public:
    Value(YAML::Node n) : node_(n) {
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
   * Read config files from the specified \p path.
   *
   * \param path The directory where the SMASH config files are located.
   */
  explicit Configuration(const boost::filesystem::path &path);

  /**
   * The default interface for SMASH to read configuration values.
   *
   * \param group
   * \param key
   *
   * \return A proxy object that converts to the correct type automatically on
   *         assignment.
   */
  Value take(std::initializer_list<const char *> keys) const;

  /**
   * Returns the object behind the requested \p key.
   *
   * You don't need to care about the return type of this function. Capture the
   * return value in a variable by using \c auto.
   * You can either call \c as<T>() with a specific type to convert the object
   * into a known value type, or you can use subsequent bracket operator calls
   * to get deeper in the configuration hierarchy.
   *
   * \return An opaque object where you can call operator[](const char *) or
   * as<T>().
   *
   * \see take
   */
  decltype(config_general_[""]) operator[](const char *key) const {
    return config_general_[key];
  }
};
#endif  // SRC_INCLUDE_CONFIGURATION_H_
