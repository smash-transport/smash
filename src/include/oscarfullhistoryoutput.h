/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_OSCARFULLHISTORYOUTPUT_H_
#define SRC_INCLUDE_OSCARFULLHISTORYOUTPUT_H_

#include "outputinterface.h"

#include "filedeleter.h"
#include "forwarddeclarations.h"
#include <string>

namespace Smash {

class OscarFullHistoryOutput : public OutputInterface {
 public:
  OscarFullHistoryOutput(bf::path path, Options op);
  ~OscarFullHistoryOutput();

  /// writes the initial particle information of an event
  void at_eventstart(const Particles &particles,
                     const int event_number) override;

  /// writes the final particle information of an event
  void at_eventend(const Particles &particles, const int event_number) override;

  /**
   * Write a prefix line and a line per particle to OSCAR output.
   */
  void write_interaction(const ParticleList &incoming_particles,
                         const ParticleList &outgoing_particles) override;

  void after_Nth_timestep(const Particles &particles, const int event_number,
                          const Clock &clock) override;

 protected:
  OscarFullHistoryOutput(bf::path path, const char* format_specifier,
                         Options op);
  void write(const Particles &particles);
  void write_2013(const Particles &particles);
  FilePtr file_;
  /* Use 2013 format */
  bool modern_format_;

 private:
  /* Output options */
  /* Print  initial and final particles */
  bool write_initial_final_lists_;
};
}  // namespace Smash

#endif  // SRC_INCLUDE_OSCARFULLHISTORYOUTPUT_H_
