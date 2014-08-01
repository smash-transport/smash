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

#include "oscaroutput.h"
#include "filedeleter.h"
#include "forwarddeclarations.h"
#include "configuration.h"

namespace Smash {

/**
 * \ingroup output
 */
class OscarFullHistoryOutput : public OscarOutput {
 public:
  OscarFullHistoryOutput(bf::path path, Configuration&& conf);
  ~OscarFullHistoryOutput();

  /// writes the initial particle information of an event
  void at_eventstart(const Particles &particles,
                     const int event_number) override;

  /// writes the final particle information of an event
  void at_eventend(const Particles &particles, const int event_number) override;

  /**
   * Write a prefix line and a line per particle to OSCAR output.
   */
  void at_interaction(const ParticleList &incoming_particles,
                         const ParticleList &outgoing_particles) override;

 private:
  /* Print  initial and final particles */
  bool print_start_end_;
};
}  // namespace Smash

#endif  // SRC_INCLUDE_OSCARFULLHISTORYOUTPUT_H_
