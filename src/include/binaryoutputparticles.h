/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_BINARYOUTPUTPARTICLES_H_
#define SRC_INCLUDE_BINARYOUTPUTPARTICLES_H_

#include <string>

#include "binaryoutputcollisions.h"
#include "filedeleter.h"
#include "forwarddeclarations.h"
#include "configuration.h"

namespace Smash {

/**
 * \brief SMASH output to binary file
 * ----------------------------------------------------------------------------
 * SMASH output to binary file is similar to OSCAR output,
 * but it is stored in a binary format. Such format is faster
 * to read and write, but may be architecture dependent.
 *
 * Binary file format is documented on the wiki in User Guide section
 **/

class BinaryOutputParticles : public BinaryOutputBase {
 public:
  BinaryOutputParticles(bf::path path, Configuration&& config);

  /// writes the initial particle information of an event
  void at_eventstart(const Particles &particles,
                     const int event_number) override;

  /// writes the final particle information of an event
  void at_eventend(const Particles &particles, const int event_number) override;

  void write_interaction(const ParticleList &incoming_particles,
                         const ParticleList &outgoing_particles) override;
  /// writes particles every time interval fixed by option OUTPUT_INTERVAL
  void after_Nth_timestep(const Particles &particles, const int event_number,
                          const Clock &clock) override;

 private:
  /// Option: print initial and final particles or not
  bool only_final_;
};


}  // namespace Smash

#endif  // SRC_INCLUDE_BINARYOUTPUTPARTICLES_H_
