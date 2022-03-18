/*
 *
 *    Copyright (c) 2022-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <boost/filesystem.hpp>

#include "configuration.h"

#ifndef SRC_INCLUDE_SMASH_LIBRARY_H_
#define SRC_INCLUDE_SMASH_LIBRARY_H_

namespace smash {

/* Free functions to interface with smash as a library,
 * also used in smash main function.*/

/// Set up configuration and logging from input files and extra config
Configuration setup_config_and_logging(
    const bf::path &config_file, const bf::path &particles_file = {},
    const bf::path &decaymodes_file = {},
    const std::vector<std::string> &extra_config = {});

/**
 * Initialize the particles and decays from the configuration, plus tabulate
 * the resonance integrals.
 */
void initalize_particles_decays_and_tabulations(Configuration &configuration,
                                                std::string version,
                                                bf::path tabulations_path = {});

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_LIBRARY_H_
