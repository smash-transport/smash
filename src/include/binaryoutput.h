/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_BINARYOUTPUT_H_
#define SRC_INCLUDE_BINARYOUTPUT_H_

#include "outputinterface.h"

#include "filedeleter.h"
#include "forwarddeclarations.h"

namespace Smash {

class BinaryOutput : public OutputInterface {
 public:
  BinaryOutput(bf::path path);
  ~BinaryOutput();

  /// writes the initial particle information of an event
  void at_eventstart(const Particles &particles,
                     const int event_number) override;

  /// writes the final particle information of an event
  void at_eventend(const Particles &particles, const int event_number) override;

  void write_interaction(const ParticleList &incoming_particles,
                         const ParticleList &outgoing_particles) override;

  void after_Nth_timestep(const Particles &particles, const int event_number,
                          const Clock &clock) override;

 protected:
  
  void write(const Particles &particles, const int event_number);
  FilePtr file_;
  const char* smash_version_string = "SMASH V0.3\n";
  const char* format_version_string = "BINARY TEST 0.1\n";
  const char* contents_units = "id PDGid t[fm/c] x[fm] y[fm] z[fm] E[GeV] px[GeV/c] py[GeV/c] pz[GeV/c] int(reserved) int(reserved) double(reserved)\n";

};
}  // namespace Smash

#endif  // SRC_INCLUDE_BINARYOUTPUT_H_
