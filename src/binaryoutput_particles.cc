/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/binaryoutput_particles.h"

#include <string>
#include <boost/filesystem.hpp>
#include <include/config.h>

#include "include/clock.h"
#include "include/forwarddeclarations.h"
#include "include/particles.h"
#include "include/inputfunctions.h"
#include "include/configuration.h"

namespace Smash {

BinaryOutputParticles::BinaryOutputParticles(bf::path path,
                                              Configuration&& config)
    : file_{std::fopen(((
            path / "particles_binary.bin")).native().c_str(), "wb")},
  only_final_(config.has_value({"only_final"}) 
                           ? config.take({"only_final"}) : true) {

  fwrite("SMSH", 4, 1, file_.get());  // magic number
  write(0);              // file format version number
  write(std::to_string(VERSION_MAJOR));  // version
}

void BinaryOutputParticles::at_eventstart(const Particles &particles,
                                 const int /*event_number*/) {
  char pchar = 'p';
  if (!only_final_) {
    std::fwrite(&pchar, sizeof(char), 1, file_.get());
    write(particles.size());
    write(particles);
  }
}

void BinaryOutputParticles::at_eventend(const Particles &particles,
                               const int event_number) {
  char pchar = 'p';
  if (only_final_) {
    std::fwrite(&pchar, sizeof(char), 1, file_.get());
    write(particles.size());
    write(particles);
  }

  // Event end line
  char fchar = 'f';
  std::fwrite(&fchar, sizeof(char), 1, file_.get());
  write(event_number);
  
  /* Flush to disk */
  std::fflush(file_.get());
}

void BinaryOutputParticles::write_interaction(const ParticleList &/*incoming*/,
                                     const ParticleList &/*outgoing*/) {
  /* No output of this kind in particles output */
}

void BinaryOutputParticles::after_Nth_timestep(const Particles &particles,
                                      const int /*event_number*/,
                                      const Clock &) {
  char pchar = 'p';
  if (!only_final_) {
    std::fwrite(&pchar, sizeof(char), 1, file_.get());
    write(particles.size());
    write(particles);
  }
}

// write functions:
inline void BinaryOutputParticles::write(const std::string &s) {
  std::int32_t size = s.size();
  std::fwrite(&size, sizeof(std::int32_t), 1, file_.get());
  std::fwrite(s.c_str(), s.size(), 1, file_.get());
}

inline void BinaryOutputParticles::write(const FourVector &v) {
  std::fwrite(v.begin(), sizeof(*v.begin()), 4, file_.get());
}

inline void BinaryOutputParticles::write(std::int32_t x) {
  std::fwrite(&x, sizeof(x), 1, file_.get());
}

void BinaryOutputParticles::write(const Particles &particles) {
  for (const auto &p : particles.data()) {
    write(p.momentum());
    write(p.position());
    write(p.pdgcode().get_decimal());
    write(p.id());
  }
}

}  // namespace Smash
