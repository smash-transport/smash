/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_OSCAROUTPUT_H_
#define SRC_INCLUDE_OSCAROUTPUT_H_

#include "filedeleter.h"
#include "forwarddeclarations.h"
#include "outputinterface.h"
#include <boost/filesystem.hpp>

namespace Smash {
class Particles;

class OscarOutput : public OutputInterface {
 public:
  OscarOutput(bf::path path);
  ~OscarOutput();

  /// writes the initial particle information of an event
  void at_eventstart(const Particles &particles, const int event_number) override;
  /// writes the final particle information of an event
  void at_eventend(const Particles &particles, const int event_number) override;

  /**
   * Write a prefix line and a line per particle to OSCAR output.
   */
  void write_interaction(const ParticleList &incoming_particles,
                         const ParticleList &outgoing_particles) override;
  void after_Nth_timestep(const Particles &particles, const int event_number, const int timestep) override;

 private:
  void write(const Particles &particles);

  FilePtr file_;
};
}  // namespace Smash

#endif  // SRC_INCLUDE_OSCAROUTPUT_H_
