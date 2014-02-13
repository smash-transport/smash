/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/configuration.h"

#include <cstdio>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>
#include <yaml-cpp/yaml.h>

namespace bf = boost::filesystem;

Configuration::Configuration(const bf::path &path) {
  const auto file_path = path / "config_general.yaml";
  assert(bf::exists(file_path));
  config_general_ = YAML::LoadFile(file_path.native());
}
