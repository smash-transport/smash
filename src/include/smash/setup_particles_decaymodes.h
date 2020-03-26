/*
 *
 *    Copyright (c) 2019-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SETUP_PARTICLES_DECAYMODES_H_
#define SRC_INCLUDE_SETUP_PARTICLES_DECAYMODES_H_

#include <string>
#include <utility>

namespace smash {
/**
 * Loads particles and decaymodes from provided files
 * particles_file and decaymodes_file. In case if particles_file
 * or decaymodes_file are nullptr, the defaults are taken
 * \param[in] particles_file a file containing particles list.
 *            See \ref inputparticles.
 * \param[in] decaymodes_file a file containing decay modes of
 *            the resonances. See \ref inputdecaymodes.
 * \return a pair of strings -- the contents of particle
 *             and decaymode files.
 */
std::pair<std::string, std::string> load_particles_and_decaymodes(
    const char *particles_file, const char *decaymodes_file);
/// Loads default smash particle list and decaymodes
void load_default_particles_and_decaymodes();

}  // namespace smash

#endif  // SRC_INCLUDE_SETUP_PARTICLES_DECAYMODES_H_
