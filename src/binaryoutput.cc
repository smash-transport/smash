/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/binaryoutput.h"

#include <include/config.h>
#include <boost/filesystem.hpp>
#include <string>

#include "include/clock.h"
#include "include/forwarddeclarations.h"
#include "include/particles.h"
#include "include/inputfunctions.h"

namespace Smash {

BinaryOutput::BinaryOutput(bf::path path, Options op)
    : particles_file_{
         std::fopen((path / "particles_binary.bin").native().c_str(), "wb")},
      collisions_file_{
         std::fopen((path / "collisions_binary.bin").native().c_str(), "wb")} {
  enable_collision_output_ = str_to_bool(op["collisions_output"]);
  only_final_ = str_to_bool(op["only_final"]);
  print_start_end_ = str_to_bool(op["print_start_end"]);

  if (enable_collision_output_) {
    fwrite("SMSH", 4, 1, collisions_file_.get());  // magic number
    write(0, collisions_file_.get());              // file format version number
    write(std::to_string(VERSION_MAJOR), collisions_file_.get());  // version

    write(10, collisions_file_.get());  // number of records
    write(
        "id PDGid t[fm/c] x[fm] y[fm] z[fm] E[GeV] px[GeV/c] py[GeV/c] "
        "pz[GeV/c]",
        collisions_file_.get());
  }

  fwrite("SMSH", 4, 1, particles_file_.get());  // magic number
  write(0, particles_file_.get());              // file format version number
  write(std::to_string(VERSION_MAJOR),
        particles_file_.get());  // version of SMASH

  // the following header information is unnecessary because it is implicitly
  // defined by the "file format version number" above already.
  write(10, particles_file_.get());  // number of records
  write(
      "id PDGid t[fm/c] x[fm] y[fm] z[fm] E[GeV] px[GeV/c] py[GeV/c] "
      "pz[GeV/c]",
      particles_file_.get());
}

void BinaryOutput::at_eventstart(const Particles &particles,
                                 const int event_number) {
  const size_t zero = 0;
  if (enable_collision_output_ && print_start_end_) {
    write(zero, collisions_file_.get());
    write(particles.size(), collisions_file_.get());
    write(particles, collisions_file_.get());
  }

  if (!only_final_) {
    write(particles.size(), particles_file_.get());
    write(event_number, particles_file_.get());
    write(particles, particles_file_.get());
  }
}

void BinaryOutput::at_eventend(const Particles &particles,
                               const int event_number) {
  const size_t zero = 0;
  if (enable_collision_output_) {
    if (print_start_end_) {
      write(particles.size(), collisions_file_.get());
      write(zero, collisions_file_.get());
      write(particles, collisions_file_.get());
    }
    // (nin, nout) = (0, 0)
    write(zero, collisions_file_.get());
    write(zero, collisions_file_.get());
    write(event_number, collisions_file_.get());
  }

  if (only_final_) {
    write(particles.size(), particles_file_.get());
    write(event_number, particles_file_.get());
    write(particles, particles_file_.get());
  }
  /* (nin, nout) = (0, 0) to indicate event end */
  write(zero, particles_file_.get());
  write(zero, particles_file_.get());
  write(event_number, particles_file_.get());

  /* Flush to disk */
  std::fflush(particles_file_.get());
  std::fflush(collisions_file_.get());
}

void BinaryOutput::write_interaction(const ParticleList &incoming,
                                     const ParticleList &outgoing) {
  if (enable_collision_output_) {
    write(incoming.size(), collisions_file_.get());
    write(outgoing.size(), collisions_file_.get());
    write(incoming, collisions_file_.get());
    write(outgoing, collisions_file_.get());
  }
}

void BinaryOutput::after_Nth_timestep(const Particles &particles,
                                      const int event_number,
                                      const Clock &) {
  if (!only_final_) {
    write(particles.size(), particles_file_.get());
    write(event_number, particles_file_.get());
    write(particles, particles_file_.get());
  }
}

// write functions:
inline void BinaryOutput::write(const std::string &s, FILE *file) {
  std::int32_t size = s.size();
  std::fwrite(&size, sizeof(std::int32_t), 1, file);
  std::fwrite(s.c_str(), s.size(), 1, file);
}

inline void BinaryOutput::write(const FourVector &v, FILE *file) {
  std::fwrite(v.begin(), sizeof(*v.begin()), 4, file);
}

inline void BinaryOutput::write(std::int32_t x, FILE *file) {
  std::fwrite(&x, sizeof(x), 1, file);
}

void BinaryOutput::write(const Particles &particles, FILE *file) {
  for (const auto &p : particles.data()) {
    write(p.id(), file);
    write(p.pdgcode().get_decimal(), file);
    write(p.position(), file);
    write(p.momentum(), file);
  }
}

void BinaryOutput::write(const ParticleList &particles,
                         FILE *file) {
  for (const auto &p : particles) {
    write(p.id(), file);
    write(p.pdgcode().get_decimal(), file);
    write(p.position(), file);
    write(p.momentum(), file);
  }
}

}  // namespace Smash
