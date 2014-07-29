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

#include <string>

#include "outputinterface.h"
#include "filedeleter.h"
#include "forwarddeclarations.h"
#include "configuration.h"

namespace Smash {

class OscarOutput : public OutputInterface {
 public:
  OscarOutput(bf::path path, std::string filename, Configuration&& conf);
  ~OscarOutput();

  void at_eventstart(const Particles &particles,
                     const int event_number) override;

  void at_eventend(const Particles &particles, const int event_number) override;

  void at_interaction(const ParticleList &incoming_particles,
                         const ParticleList &outgoing_particles) override;

  void at_intermediate_time(const Particles &particles, const int event_number,
                          const Clock &clock) override;

 protected:
  void write_format_description(void);
  void write_particledata(const ParticleData &data);
  void write(const Particles &particles);
  FilePtr file_;
  bool modern_format_;
};
}  // namespace Smash

#endif  // SRC_INCLUDE_OSCAROUTPUT_H_
