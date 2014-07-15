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
    fwrite("SMSH", 4, 1, collisions_file_.get());            // magic number
    write(0, "collisions");           // file format version number
    write(std::to_string(VERSION_MAJOR), "collisions");  // version

    write(10, "collisions");  // number of records
    write(
       "id PDGid t[fm/c] x[fm] y[fm] z[fm] E[GeV] px[GeV/c] py[GeV/c] "
       "pz[GeV/c]", "collisions");
  }

  fwrite("SMSH", 4, 1, particles_file_.get());            // magic number
  write(0, "particles");                 // file format version number
  write(std::to_string(VERSION_MAJOR), "particles");  // version of SMASH

  // the following header information is unnecessary because it is implicitly
  // defined by the "file format version number" above already.
  write(10, "particles");  // number of records
  write(
      "id PDGid t[fm/c] x[fm] y[fm] z[fm] E[GeV] px[GeV/c] py[GeV/c] "
      "pz[GeV/c]", "particles");
}

void BinaryOutput::at_eventstart(const Particles &particles,
                                 const int event_number) {
  const size_t zero = 0;
  if (enable_collision_output_ && print_start_end_) {
    write(zero, "collisions");
    write(particles.size(), "collisions");
    write(particles, "collisions");
  }

  if (!only_final_) {
    write(particles.size(), "particles");
    write(event_number, "particles");
    write(particles, "particles");
  }
}

void BinaryOutput::at_eventend(const Particles &particles,
                               const int event_number) {
  const size_t zero = 0;
  if (enable_collision_output_) {
    if (print_start_end_) {
      write(particles.size(), "collisions");
      write(zero, "collisions");
      write(particles, "collisions");
    }
    write(zero, "collisions");
    write(zero, "collisions");
    write(event_number, "collisions");
  }

  if (only_final_) {
    write(particles.size(), "particles");
    write(event_number, "particles");
    write(particles, "particles");
  }
  /* Flush to disk */
  std::fflush(particles_file_.get());
  std::fflush(collisions_file_.get());
}

void BinaryOutput::write_interaction(const ParticleList &incoming,
                                     const ParticleList &outgoing) {
  if (enable_collision_output_) {
    write(incoming.size(), "collisions");
    write(outgoing.size(), "collisions");
    write(incoming, "collisions");
    write(outgoing, "collisions");
  }
}

void BinaryOutput::after_Nth_timestep(const Particles &particles,
                                      const int event_number,
                                      const Clock &) {
  if (!only_final_) {
    write(particles.size(), "particles");
    write(event_number, "particles");
    write(particles, "particles");
  }
}

// write functions:
inline void BinaryOutput::write(const std::string &s, const std::string &op) {
  std::int32_t size = s.size();
  if (op == "collisions") {
    std::fwrite(&size, sizeof(std::int32_t), 1, collisions_file_.get());
    std::fwrite(s.c_str(), s.size(), 1, collisions_file_.get());
  } else {
    std::fwrite(&size, sizeof(std::int32_t), 1, particles_file_.get());
    std::fwrite(s.c_str(), s.size(), 1, particles_file_.get());
  }
}

inline void BinaryOutput::write(const FourVector &v, const std::string &op) {
  if (op == "collisions") {
    std::fwrite(v.begin(), sizeof(*v.begin()), 4, collisions_file_.get());
  } else {
    std::fwrite(v.begin(), sizeof(*v.begin()), 4, particles_file_.get());
  }
}

inline void BinaryOutput::write(std::int32_t x, const std::string &op) {
  if (op == "collisions") {
    std::fwrite(&x, sizeof(x), 1, collisions_file_.get());
  } else {
    std::fwrite(&x, sizeof(x), 1, particles_file_.get());
  }
}

void BinaryOutput::write(const Particles &particles, const std::string &op) {
  for (const auto &p : particles.data()) {
    write(p.id(), op);
    write(p.pdgcode().get_decimal(), op);
    write(p.position(), op);
    write(p.momentum(), op);
  }
}

void BinaryOutput::write(const ParticleList &particles,
                         const std::string &op) {
  for (const auto &p : particles) {
    write(p.id(), op);
    write(p.pdgcode().get_decimal(), op);
    write(p.position(), op);
    write(p.momentum(), op);
  }
}

}  // namespace Smash
