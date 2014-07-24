/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/binaryoutputcollisions.h"

#include <string>
#include <boost/filesystem.hpp>
#include <include/config.h>

#include "include/clock.h"
#include "include/forwarddeclarations.h"
#include "include/particles.h"
#include "include/inputfunctions.h"
#include "include/configuration.h"

namespace Smash {

BinaryOutputCollisions::BinaryOutputCollisions(bf::path path,
                                               Configuration &&config)
    : BinaryOutputBase(std::fopen(
          ((path / "collisions_binary.bin")).native().c_str(), "wb")),
      print_start_end_(config.has_value({"print_start_end"})
                           ? config.take({"print_start_end"})
                           : false) {
  fwrite("SMSH", 4, 1, file_.get());  // magic number
  write(0);              // file format version number
  write(VERSION_MAJOR);  // SMASH version
}

void BinaryOutputCollisions::at_eventstart(const Particles &particles,
                                 const int /*event_number*/) {
  char pchar = 'p';
  if (print_start_end_) {
    std::fwrite(&pchar, sizeof(char), 1, file_.get());
    write(particles.size());
    write(particles);
  }
}

void BinaryOutputCollisions::at_eventend(const Particles &particles,
                               const int event_number) {
  char pchar = 'p';
  if (print_start_end_) {
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

void BinaryOutputCollisions::write_interaction(const ParticleList &incoming,
                                     const ParticleList &outgoing) {
  char ichar = 'i';
  std::fwrite(&ichar, sizeof(char), 1, file_.get());
  write(incoming.size());
  write(outgoing.size());
  write(incoming);
  write(outgoing);
}

void BinaryOutputCollisions::after_Nth_timestep(const Particles &/*particles*/,
                                      const int /*event_number*/,
                                      const Clock &) {
  /* No output of this kind in collisions output */
}

// write functions:
void BinaryOutputBase::write(const std::string &s) {
  std::int32_t size = s.size();
  std::fwrite(&size, sizeof(std::int32_t), 1, file_.get());
  std::fwrite(s.c_str(), s.size(), 1, file_.get());
}

void BinaryOutputBase::write(const FourVector &v) {
  std::fwrite(v.begin(), sizeof(*v.begin()), 4, file_.get());
}

void BinaryOutputBase::write(const Particles &particles) {
  for (const auto &p : particles.data()) {
    write(p.momentum());
    write(p.position());
    write(p.pdgcode().get_decimal());
    write(p.id());
  }
}

void BinaryOutputBase::write(const ParticleList &particles) {
  for (const auto &p : particles) {
    write(p.momentum());
    write(p.position());
    write(p.pdgcode().get_decimal());
    write(p.id());
  }
}

}  // namespace Smash
