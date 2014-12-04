/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_DENSITYOUTPUT_H_
#define SRC_INCLUDE_DENSITYOUTPUT_H_

#include <string>

#include "filedeleter.h"
#include "outputinterface.h"
#include "forwarddeclarations.h"
#include "configuration.h"
#include "threevector.h"

namespace Smash {

/**
 * \ingroup output
 *
 * \brief Writes the density at a specified point
 *
 **/
class DensityOutput : public OutputInterface {
 public:

  DensityOutput(bf::path path, Configuration&& conf,
                const double sigma, const int ntest);
  ~DensityOutput();

  /// writes the initial particle information of an event
  void at_eventstart(const Particles &particles,
                     const int event_number) override;

  /// writes the final particle information of an event
  void at_eventend(const Particles &particles, const int event_number) override;

  /// writes particles every time interval fixed by option Output_Interval
  void at_intermediate_time(const Particles &particles, const int event_number,
                          const Clock &clock) override;
 private:
  FilePtr file_;

  /// Point, where density is calculated
  const ThreeVector r_;

  /// Sigma of Gaussian smearing for density
  const double sigma_;

  /// Number of testparticles
  const int ntest_;
};


}  // namespace Smash

#endif  // SRC_INCLUDE_DENSITYOUTPUT_H_
