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

#include <string>

#include <yaml-cpp/yaml.h>

namespace boost {
namespace filesystem {
class path;
}  // namespace filesystem
}  // namespace boost

class Configuration {
  YAML::Node config_general_;

 public:
  explicit Configuration(const boost::filesystem::path &path);

  decltype(config_general_[""]) operator[](const char *key) {
    return config_general_[key];
  }
};
#endif  // SRC_INCLUDE_CONFIGURATION_H_
