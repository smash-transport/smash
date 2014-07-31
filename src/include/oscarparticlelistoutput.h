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

#include "filedeleter.h"
#include "forwarddeclarations.h"
#include "configuration.h"
#include "oscaroutput.h"

namespace Smash {

class OscarParticleListOutput : public OscarOutput {
 public:
  OscarParticleListOutput(bf::path path, Configuration&& conf);
  ~OscarParticleListOutput();

  /// writes the initial particle information of an event
  void at_eventstart(const Particles &particles,
                     const int event_number) override;

  /// writes the final particle information of an event
  void at_eventend(const Particles &particles, const int event_number) override;

  void at_intermediate_time(const Particles &particles, const int event_number,
                          const Clock &clock) override;

 private:
  /// An option. If true - only final particles in event are printed
  bool only_final_;
};
}  // namespace Smash

#endif  // SRC_INCLUDE_OSCARPARTICLELISTOUTPUT_H_
