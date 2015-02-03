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
 * \ingroup output
 *
 * \brief Writes the particle list at specific times to binary file
 *
 * This class writes the current particle list at a specific time t
 * to the binary output file. These specific time can
 * be: event start, event end, every next time interval \f$\Delta t \f$.
 * Writing (or not writing) output at these moments is controlled by options.
 * Time interval \f$\Delta t \f$ is also regulated by an option.
 * Output file is binary and has a block structure.
 *
 * Details of the output format can be found
 * on the wiki in User Guide section, look for binary output.
 **/
class BinaryOutputParticles : public BinaryOutputBase {
 public:
  BinaryOutputParticles(bf::path path, Configuration&& config);

  /// writes the initial particle information of an event
  void at_eventstart(const Particles &particles,
                     const int event_number) override;

  /// writes the final particle information of an event
  void at_eventend(const Particles &particles, const int event_number) override;

  void at_interaction(const ParticleList &incoming_particles,
                      const ParticleList &outgoing_particles,
                      const double density,
                      const double total_cross_section) override;
  /// writes particles every time interval fixed by option OUTPUT_INTERVAL
  void at_intermediate_time(const Particles &particles, const int event_number,
                          const Clock &clock) override;

 private:
  /// Option: print initial and final particles or not
  bool only_final_;
};


}  // namespace Smash

#endif  // SRC_INCLUDE_BINARYOUTPUTPARTICLES_H_
