/*
 *
 *    Copyright (c) 2022,2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <string>
#include <vector>

#include "configuration.h"
#include "sha256.h"

#ifndef SRC_INCLUDE_SMASH_LIBRARY_H_
#define SRC_INCLUDE_SMASH_LIBRARY_H_

namespace smash {

/* Free functions to interface with smash as a library,
 * also used in smash main function. */

/**
 * Set up configuration and logging from input files and extra config
 *
 * \param[in] config_file Path to config input file
 * \param[in] particles_file Path to particles input file.
 * \param[in] decaymodes_file Path to  decaymodes input file.
 * \param[in] extra_config Extra config entries.
 *
 * \return Configuration object with particles, decaymodes and extra
 * configs included.
 *
 * If no particles and decaymodes files are given the default files in
 * the input directory are used.
 */
Configuration setup_config_and_logging(
    const std::string &config_file, const std::string &particles_file = {},
    const std::string &decaymodes_file = {},
    const std::vector<std::string> &extra_config = {});

/**
 * Wrapper over a function that initializes the particles and decays from the
 * given configuration, and over another that tabulates the resonance integrals.
 *
 * \param[in] configuration Fully-setup configuration i.e. including
 * particles and decaymodes.
 * \param[in] version Current version of SMASH.
 * \param[in] tabulations_dir Path where tabulations should be stored.
 */
void initialize_particles_decays_and_tabulations(
    Configuration &configuration, const std::string &version,
    const std::string &tabulations_dir = {});

/**
 * Initialize the particles and decays from the given configuration.
 *
 * \param[in] configuration Fully-setup configuration i.e. including
 * particles and decaymodes.
 * \param[in] version Current version of SMASH.
 * \return hash of the version, particle list, and decay modes.
 */
sha256::Hash initialize_particles_decays_and_return_hash(
    Configuration &configuration, const std::string &version);

/**
 * Tabulate the resonance integrals.
 *
 * \param[in] hash Hash of the SMASH version, particle list, and decay modes.
 * \param[in] tabulations_dir Path where tabulations should be stored.
 */
void tabulate_resonance_integrals(const sha256::Hash &hash,
                                  const std::string &tabulations_dir);
}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_LIBRARY_H_
