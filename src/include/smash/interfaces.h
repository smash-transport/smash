
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


#ifndef SRC_INCLUDE_SMASH_SMASHINTERFACES_H_
#define SRC_INCLUDE_SMASH_SMASHINTERFACES_H_

namespace smash {

/// Set up configuration and logging from input files and extra config
Configuration configure(const bf::path &config_file,
                        const char *particles_file = nullptr,
                        const char *decaymodes_file = nullptr,
                        const std::vector<std::string> &extra_config = {});



/**
 * Initialize the particles and decays from the configuration, plus tabulate
 * the resonance integrals.
 */
void initalize(Configuration &configuration, std::string version,
               bf::path tabulations_path);




////////////////////////////////////////////////////////////////////////////////





Configuration configure(const bf::path &config_file,
                        const char *particles_file,
                        const char *decaymodes_file ,
                        const std::vector<std::string> &extra_config) {
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
  if (configuration.has_value({"particles"}) && particles_file) {
    logg[LMain].warn(
        "Ambiguity: particles from external file ", particles_file,
        " requested, but there is also particle list in the config."
        " Using particles from ",
        particles_file);
  }
  if (!configuration.has_value({"particles"}) || particles_file) {
    configuration["particles"] = particles_and_decays.first;
  }

  if (configuration.has_value({"decaymodes"}) && decaymodes_file) {
    logg[LMain].warn(
        "Ambiguity: decaymodes from external file ", decaymodes_file,
        " requested, but there is also decaymodes list in the config."
        " Using decaymodes from",
        decaymodes_file);
  }
  if (!configuration.has_value({"decaymodes"}) || decaymodes_file) {
    configuration["decaymodes"] = particles_and_decays.second;
  }

  // // TODO(stdnmr) Something fishy might go on here, need to investigate
  // santizer tests say: smash-devel/src/smash.cc:372:15: runtime error: execution reached the end of a value-returning function without returning a value
  // SUMMARY: UndefinedBehaviorSanitizer: undefined-behavior /Users/stdnmr/smash-devel/src/smash.cc:372:15 in
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

#endif  // SRC_INCLUDE_SMASH_SMASHINTERFACES_H_
