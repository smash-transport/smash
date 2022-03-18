/*
 *
 *    Copyright (c) 2022-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <string>
#include <vector>

#include "configuration.h"

#ifndef SRC_INCLUDE_SMASH_LIBRARY_H_
#define SRC_INCLUDE_SMASH_LIBRARY_H_

namespace smash {

/* Free functions to interface with smash as a library,
 * also used in smash main function.*/

/// Set up configuration and logging from input files and extra config
Configuration setup_config_and_logging(
    const std::string &config_file, const std::string &particles_file = {},
    const std::string &decaymodes_file = {},
    const std::vector<std::string> &extra_config = {});

/**
 * Initialize the particles and decays from the configuration, plus tabulate
 * the resonance integrals.
 */
void initalize_particles_decays_and_tabulations(
    Configuration &configuration, const std::string &version,
    const std::string &tabulations_dir = {});

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_LIBRARY_H_
