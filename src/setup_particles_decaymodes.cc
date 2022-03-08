/*
 *
 *    Copyright (c) 2019-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/setup_particles_decaymodes.h"

#include "smash/decaymodes.h"
#include "smash/inputfunctions.h"

namespace {
#ifndef DOXYGEN
namespace particles_txt {
#include <particles.txt.h>
}  // namespace particles_txt
namespace decaymodes_txt {
#include <decaymodes.txt.h>
}  // namespace decaymodes_txt
#endif
}  // unnamed namespace

namespace smash {

std::pair<std::string, std::string> load_particles_and_decaymodes(
    const boost::filesystem::path &particles_file,
    const boost::filesystem::path &decaymodes_file) {
  std::string particle_string, decay_string;
  if (!particles_file.empty()) {
    if (!boost::filesystem::exists(particles_file)) {
      std::stringstream err;
      err << "The particles file was expected at '" << particles_file
          << "', but the file does not exist.";
      throw std::runtime_error(err.str());
    }
    particle_string = read_all(boost::filesystem::ifstream{particles_file});
    if (has_crlf_line_ending(particle_string)) {
      std::stringstream err;
      err << "The particles file has CR LF line endings. Please use LF"
             " line endings.";
      throw std::runtime_error(err.str());
    }
  } else {
    particle_string = particles_txt::data;
  }

  if (!decaymodes_file.empty()) {
    if (!boost::filesystem::exists(decaymodes_file)) {
      std::stringstream err;
      err << "The decay modes file was expected at '" << decaymodes_file
          << "', but the file does not exist.";
      throw std::runtime_error(err.str());
    }
    decay_string = read_all(boost::filesystem::ifstream{decaymodes_file});
    if (has_crlf_line_ending(decay_string)) {
      std::stringstream err;
      err << "The decay mode file has CR LF line endings. Please use LF"
             " line endings.";
      throw std::runtime_error(err.str());
    }
  } else {
    decay_string = decaymodes_txt::data;
  }
  return std::make_pair(particle_string, decay_string);
}

void intialize_default_particles_and_decaymodes() {
  const auto pd = load_particles_and_decaymodes({}, {});
  ParticleType::create_type_list(pd.first);
  DecayModes::load_decaymodes(pd.second);
  ParticleType::check_consistency();
}

}  // namespace smash
