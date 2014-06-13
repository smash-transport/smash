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

#include <boost/filesystem.hpp>

namespace Smash {

BinaryOutput::BinaryOutput(bf::path path)
  : file_{std::fopen((path / "particles_binary.bin").native().c_str(), "wb")} {

  const char smash_version_string[80] = "SMASH V0.3\n";
  const char format_version_string[80] = "BINARY TEST 0.1\n";
  const char contents_units[300] = "/2 ints:/ id PDGid /8 doubles:/"
                                   "t[fm/c] x[fm] y[fm] z[fm] E[GeV] px[GeV/c]"
                                   "py[GeV/c] pz[GeV/c] int(reserved)"
                                   "int(reserved) double(reserved)\n";
  
  char sizeof_int_string[30];
  char sizeof_double_string[30];
   
  snprintf(sizeof_int_string, sizeof(sizeof_int_string),
           "size of int is %zu bytes\n", sizeof(int));
  snprintf(sizeof_double_string, sizeof(sizeof_double_string),
           "size of double is %zu bytes\n", sizeof(double));

  fwrite (&smash_version_string , sizeof(char), 80, file_.get());
  fwrite (&format_version_string , sizeof(char), 80, file_.get());
  fwrite (&contents_units, sizeof(char), 300, file_.get());
  fwrite (&sizeof_int_string, sizeof(char), 30, file_.get());
  fwrite (&sizeof_double_string, sizeof(char), 30, file_.get());
}


BinaryOutput::~BinaryOutput() {}

void BinaryOutput::at_eventstart(const Particles &/*particles*/,
                                const int /*event_number*/) { 
}

void BinaryOutput::at_eventend(const Particles &particles,
                              const int event_number) {
 
  write(particles, event_number);

  /* Flush to disk */
  std::fflush(file_.get());
}

void BinaryOutput::write(const Particles &particles, const int event_number) {

  double r[4], mom[4], dummy_double;
  int id, pdgcode, dummy_int[2];
  int psize = static_cast<int>(particles.size());

  fwrite (&psize,         sizeof(int), 1, file_.get());
  fwrite (&event_number,  sizeof(int), 1, file_.get());

  for (const auto &p : particles.data()) {

    r[0] = p.position().x0();
    r[1] = p.position().x1();
    r[2] = p.position().x2();
    r[3] = p.position().x3();

    mom[0] = p.momentum().x0();
    mom[1] = p.momentum().x1();
    mom[2] = p.momentum().x2();
    mom[3] = p.momentum().x3();

    pdgcode = p.pdgcode().get_decimal();
    id = p.id();
 
    fwrite (&id,       sizeof(int), 1, file_.get());
    fwrite (&pdgcode,  sizeof(int), 1, file_.get());

    fwrite (r,   sizeof(double), 4, file_.get());
    fwrite (mom, sizeof(double), 4, file_.get());

    // Reserved fields
    fwrite (dummy_int, sizeof(int), 2, file_.get());
    fwrite (&dummy_double, sizeof(double), 1, file_.get());
  }

}

void BinaryOutput::write_interaction(
  const ParticleList &/*incoming_particles*/,
  const ParticleList &/*outgoing_particles*/) {
}

void BinaryOutput::after_Nth_timestep(const Particles & /*particles*/,
                                                const int /*event_number*/,
                                     const Clock& /*clock*/) {
}

}  // namespace Smash
