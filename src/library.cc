/*
 *
 *    Copyright (c) 2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/library.h"

#include <filesystem>

#include "smash/decaymodes.h"
#include "smash/inputfunctions.h"
#include "smash/isoparticletype.h"
#include "smash/logging.h"
#include "smash/setup_particles_decaymodes.h"
#include "smash/stringfunctions.h"

namespace smash {
static constexpr int LMain = LogArea::Main::id;

static void do_minimal_loggers_setup_for_config_validation() {
  const std::string conf_tag = LogArea::Configuration::textual();
  const std::string main_tag = LogArea::Main::textual();
  const auto size =
      conf_tag.size() > main_tag.size() ? conf_tag.size() : main_tag.size();
  logg[LogArea::Configuration::id].setAreaName(utf8::fill_both(conf_tag, size));
  logg[LogArea::Main::id].setAreaName(utf8::fill_both(main_tag, size));
}

Configuration setup_config_and_logging(
    const std::string &config_file, const std::string &particles_file,
    const std::string &decaymodes_file,
    const std::vector<std::string> &extra_config) {
  // Read in config file
  std::filesystem::path config_path(config_file);
  Configuration configuration(config_path.parent_path(),
                              config_path.filename());

  // Merge config passed via command line
  for (const auto &config : extra_config) {
    configuration.merge_yaml(config);
  }

  // Fully validate the configuration
  do_minimal_loggers_setup_for_config_validation();
  if (configuration.validate() == Configuration::Is::Invalid) {
    throw std::runtime_error("Validation of SMASH input failed.");
  }

  // Set up logging
  set_default_loglevel(
      configuration.take({"Logging", "default"}, einhard::ALL));
  create_all_loggers(configuration.extract_sub_configuration(
      {"Logging"}, Configuration::GetEmpty::Yes));

  logg[LMain].trace(SMASH_SOURCE_LOCATION, " load ParticleType and DecayModes");

  std::filesystem::path particles_path(particles_file);
  std::filesystem::path decaymodes_path(decaymodes_file);
  auto particles_and_decays =
      load_particles_and_decaymodes(particles_path, decaymodes_path);
  /* For particles and decaymodes: external file is superior to config.
   * However, warn in case of conflict.
   */
  if (configuration.has_value({"particles"}) && !particles_path.empty()) {
    logg[LMain].warn(
        "Ambiguity: particles from external file ", particles_path,
        " requested, but there is also particle list in the config."
        " Using particles from ",
        particles_path);
  }
  if (!configuration.has_value({"particles"}) || !particles_path.empty()) {
    configuration.set_value({"particles"}, particles_and_decays.first);
  }

  if (configuration.has_value({"decaymodes"}) && !decaymodes_path.empty()) {
    logg[LMain].warn(
        "Ambiguity: decaymodes from external file ", decaymodes_path,
        " requested, but there is also decaymodes list in the config."
        " Using decaymodes from",
        decaymodes_path);
  }
  if (!configuration.has_value({"decaymodes"}) || !decaymodes_path.empty()) {
    configuration.set_value({"decaymodes"}, particles_and_decays.second);
  }

  return configuration;
}

void initialize_particles_decays_and_tabulations(
    Configuration &configuration, const std::string &version,
    const std::string &tabulations_dir) {
  logg[LMain].trace(SMASH_SOURCE_LOCATION,
                    " create ParticleType and DecayModes");
  const std::string particles_string = configuration.take({"particles"});
  const std::string decaymodes_string = configuration.take({"decaymodes"});
  ParticleType::create_type_list(particles_string);
  DecayModes::load_decaymodes(decaymodes_string);
  ParticleType::check_consistency();

  // Calculate a hash of the SMASH version, the particles and decaymodes.
  sha256::Context hash_context;
  hash_context.update(version);
  hash_context.update(particles_string);
  hash_context.update(decaymodes_string);
  const auto hash = hash_context.finalize();
  logg[LMain].info() << "Config hash: " << sha256::hash_to_string(hash);

  logg[LMain].info("Tabulating cross section integrals...");
  std::filesystem::path tabulations_path(tabulations_dir);
  if (!tabulations_path.empty()) {
    // Store tabulations on disk
    std::filesystem::create_directories(tabulations_path);
    logg[LMain].info() << "Tabulations path: " << tabulations_path;
  }
  IsoParticleType::tabulate_integrals(hash, tabulations_path);
}

}  // namespace smash
