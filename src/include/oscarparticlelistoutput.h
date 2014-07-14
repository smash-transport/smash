/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_OSCARPARTICLELISTOUTPUT_H_
#define SRC_INCLUDE_OSCARPARTICLELISTOUTPUT_H_

#include "oscarfullhistoryoutput.h"

namespace Smash {

class OscarParticleListOutput : public OscarFullHistoryOutput {
 public:
  OscarParticleListOutput(bf::path path);
  ~OscarParticleListOutput();

  /// writes the initial particle information of an event
  void at_eventstart(const Particles &particles,
                     const int event_number) override;

  /**
   * Write a prefix line and a line per particle to OSCAR output.
   */
  void write_interaction(const ParticleList &incoming_particles,
                         const ParticleList &outgoing_particles) override;

};
}  // namespace Smash

#endif  // SRC_INCLUDE_OSCARPARTICLELISTOUTPUT_H_
