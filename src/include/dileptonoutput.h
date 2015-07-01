/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_DILEPTONOUTPUT_H_
#define SRC_INCLUDE_DILEPTONOUTPUT_H_

#include <string>

#include <boost/filesystem.hpp>

#include "configuration.h"
#include "filedeleter.h"
#include "forwarddeclarations.h"
#include "outputinterface.h"

namespace Smash {


class DileptonOutput : public OutputInterface {
public:
  DileptonOutput (bf::path path);

  void at_eventstart(const Particles &particles,
                     const int event_number) override;

  void at_eventend(const Particles &particles, const int event_number) override;

  void at_interaction(const ParticleList &incoming_particles,
                      const ParticleList &outgoing_particles,
                      const double density,
                      const double total_cross_section,
                      const ProcessType process_type) override {};

  void at_intermediate_time(const Particles &, const int,
                            const Clock &) override {};

  void dileptons(const ParticleList &incoming_particles,
                 const ParticleList &outgoing_particles,
                 float shining_weight) override;

private:
  FilePtr file_;

};

}  // namespace Smash

#endif  // SRC_INCLUDE_DILEPTONOUTPUT_H_
