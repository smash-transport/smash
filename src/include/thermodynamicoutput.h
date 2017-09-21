/*
 *
 *    Copyright (c) 2014-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_THERMODYNAMICOUTPUT_H_
#define SRC_INCLUDE_THERMODYNAMICOUTPUT_H_

#include <set>
#include <string>

#include "configuration.h"
#include "density.h"
#include "experimentparameters.h"
#include "filedeleter.h"
#include "forwarddeclarations.h"
#include "outputinterface.h"
#include "threevector.h"

namespace Smash {

/**
 * \ingroup output
 *
 * \brief Writes the thermodynamic quantities at a specified point versus time
 *
 * This class is a temporary solution to write thermodynamic
 * quantities out. Calculations are called directly inside the
 * output functions. In future it should be substituted by some
 * more general output.
 *
 **/
class ThermodynamicOutput : public OutputInterface {
 public:
  ThermodynamicOutput(const bf::path &path, Configuration &&conf);
  ~ThermodynamicOutput();

  /// writes the initial particle information of an event
  void at_eventstart(const Particles &particles,
                     const int event_number) override;

  /// writes the final particle information of an event
  void at_eventend(const Particles &particles, const int event_number) override;

  /// writes thermodynamics every time interval fixed by option Output_Interval
  void at_intermediate_time(const Particles &particles, const Clock &clock,
                            const DensityParameters &dens_param) override;

  /** Prints density along the specified line. Useful to make 1D plots of
    * density profiles.
   */
  void density_along_line(const char *file_name, const ParticleList &plist,
                          const DensityParameters &param, DensityType dens_type,
                          const ThreeVector &line_start,
                          const ThreeVector &line_end, int n_points);

 private:
  FilePtr file_;
  /// Set of quantities to be computed
  const std::set<ThermodynamicQuantity> td_set_;
  /// Point, where thermodynamic quantities are calculated
  ThreeVector r_;
  /// Type (e.g., baryon/pion/hadron) of thermodynamic quantity
  const DensityType dens_type_;
  /** Whether smearing is on or off; WARNING : if smearing is off,
      then final result is in GeV instead of GeV/fm3 */
  const bool smearing_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_THERMODYNAMICOUTPUT_H_
