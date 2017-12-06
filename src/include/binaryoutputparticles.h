/*
 *
 *    Copyright (c) 2014-2017
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
#include "outputparameters.h"

namespace smash {

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
  BinaryOutputParticles(const bf::path &path, std::string name,
                        const OutputParameters &out_par);

  /// writes the initial particle information of an event
  void at_eventstart(const Particles &particles,
                     const int event_number) override;

  /// writes the final particle information of an event
  void at_eventend(const Particles &particles, const int event_number,
                   double impact_parameter) override;

  /// writes particles every time interval fixed by option OUTPUT_INTERVAL
  void at_intermediate_time(const Particles &particles, const Clock &clock,
                            const DensityParameters &dens_param) override;

 private:
  /// Option: print initial and final particles or not
  bool only_final_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_BINARYOUTPUTPARTICLES_H_
