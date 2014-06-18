/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/binaryoutput.h"
#include "include/clock.h"
#include "include/forwarddeclarations.h"
#include "include/particles.h"
#include <include/config.h>

#include <boost/filesystem.hpp>
#include <string>

namespace Smash {

BinaryOutput::BinaryOutput(bf::path path)
    : file_{
          std::fopen((path / "particles_binary.bin").native().c_str(), "wb")} {
  fwrite("SMSH", 4, 1, file_.get());     // magic number
  write(0);                              // file format version number
  write(std::to_string(VERSION_MAJOR));  // version of SMASH

  // the following header information is unnecessary because it is implicitly
  // defined by the "file format version number" above already.
  write(10);  // number of records
  write(
      "id PDGid t[fm/c] x[fm] y[fm] z[fm] E[GeV] px[GeV/c] py[GeV/c] "
      "pz[GeV/c]");
}

void BinaryOutput::at_eventstart(const Particles &particles,
                                 const int event_number) {
  write(particles, event_number);
}

void BinaryOutput::at_eventend(const Particles &,
                               const int) {
  /* Flush to disk */
  std::fflush(file_.get());
}

void BinaryOutput::write_interaction(const ParticleList &,
                                     const ParticleList &) {}

void BinaryOutput::after_Nth_timestep(const Particles &particles,
                                      const int event_number,
                                      const Clock &) {
  write(particles, event_number);
}

// write functions:
inline void BinaryOutput::write(const std::string &s) {
  std::int32_t size = s.size();
  std::fwrite(&size, sizeof(std::int32_t), 1, file_.get());
  std::fwrite(s.c_str(), s.size(), 1, file_.get());
}

inline void BinaryOutput::write(const FourVector &v) {
  std::fwrite(v.begin(), sizeof(*v.begin()), 4, file_.get());
}

inline void BinaryOutput::write(std::int32_t x) {
  std::fwrite(&x, sizeof(x), 1, file_.get());
}

void BinaryOutput::write(const Particles &particles, const int event_number) {
  write(particles.size());
  write(event_number);

  for (const auto &p : particles.data()) {
    write(p.id());
    write(p.pdgcode().get_decimal());
    write(p.position());
    write(p.momentum());
  }
}

}  // namespace Smash
