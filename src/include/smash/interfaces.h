
/*
 *
 *    Copyright (c) 2022-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/setup_particles_decaymodes.h"
#include "smash/decaymodes.h"
#include "smash/filelock.h"
#include "smash/random.h"
#include "smash/sha256.h"
#include "smash/config.h"
#include "smash/experiment.h"
#include "smash/stringfunctions.h"

// TODO(stdnmr) check what includes are really necessary


#ifndef SRC_INCLUDE_SMASH_LIB_H_
#define SRC_INCLUDE_SMASH_LIB_H_

namespace smash {

// Free functions to interface with smash
// TODO(stdnmr) Put example program here
// TODO(stdnmr) Place implementation into .cc

/// Set up configuration and logging from input files and extra config
Configuration setup_config_and_logging(const bf::path &config_file,
                                       const bf::path &particles_file = {},
                                       const bf::path &decaymodes_file = {},
                                       const std::vector<std::string> &extra_config = {});


/**
 * Initialize the particles and decays from the configuration, plus tabulate
 * the resonance integrals.
 */
void initalize(Configuration &configuration, std::string version,
               bf::path tabulations_path);

////////////////////////////////////////////////////////////////////////////////

Configuration setup_config_and_logging(const bf::path &config_file,
                           const bf::path &particles_file,
                           const bf::path &decaymodes_file,
                           const std::vector<std::string> &extra_config){
  // Read in config file
  Configuration configuration(config_file.parent_path(),
                              config_file.filename());

  // Merge config passed via command line
  for (const auto &config : extra_config) {
    configuration.merge_yaml(config);
  }

  // Set up logging
  set_default_loglevel(
      configuration.take({"Logging", "default"}, einhard::ALL));
  create_all_loggers(configuration["Logging"]);

  logg[LMain].trace(SMASH_SOURCE_LOCATION, " load ParticleType and DecayModes");

  auto particles_and_decays =
      load_particles_and_decaymodes(particles_file, decaymodes_file);
  /* For particles and decaymodes: external file is superior to config.
   * Hovever, warn in case of conflict.
   */
  if (configuration.has_value({"particles"}) && !particles_file.empty()) {
    logg[LMain].warn(
        "Ambiguity: particles from external file ", particles_file,
        " requested, but there is also particle list in the config."
        " Using particles from ",
        particles_file);
  }
  if (!configuration.has_value({"particles"}) || !particles_file.empty()) {
    configuration["particles"] = particles_and_decays.first;
  }

  if (configuration.has_value({"decaymodes"}) && !decaymodes_file.empty()) {
    logg[LMain].warn(
        "Ambiguity: decaymodes from external file ", decaymodes_file,
        " requested, but there is also decaymodes list in the config."
        " Using decaymodes from",
        decaymodes_file);
  }
  if (!configuration.has_value({"decaymodes"}) || !decaymodes_file.empty()) {
    configuration["decaymodes"] = particles_and_decays.second;
  }

  return configuration;
}

void initalize(Configuration &configuration, std::string version,
               bf::path tabulations_path) {
  logg[LMain].trace(SMASH_SOURCE_LOCATION,
                    " create ParticleType and DecayModes");
  ParticleType::create_type_list(configuration.take({"particles"}));
  DecayModes::load_decaymodes(configuration.take({"decaymodes"}));
  ParticleType::check_consistency();

  logg[LMain].info("Tabulating cross section integrals...");
  // Calculate a hash of the SMASH version, the particles and decaymodes.
  const std::string particle_string = configuration["particles"].to_string();
  const std::string decay_string = configuration["decaymodes"].to_string();
  sha256::Context hash_context;
  hash_context.update(version);
  hash_context.update(particle_string);
  hash_context.update(decay_string);
  const auto hash = hash_context.finalize();
  logg[LMain].info() << "Config hash: " << sha256::hash_to_string(hash);
  IsoParticleType::tabulate_integrals(hash, tabulations_path);
}


} // namespace smash

#endif  // SRC_INCLUDE_SMASH_LIB_H_
