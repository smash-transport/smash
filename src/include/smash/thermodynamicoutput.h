/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_THERMODYNAMICOUTPUT_H_
#define SRC_INCLUDE_THERMODYNAMICOUTPUT_H_

#include <set>
#include <string>

#include "density.h"
#include "experimentparameters.h"
#include "file.h"
#include "forwarddeclarations.h"
#include "outputinterface.h"
#include "outputparameters.h"
#include "threevector.h"

namespace smash {

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
  ThermodynamicOutput(const bf::path &path, const std::string &name,
                      const OutputParameters &out_par);
  ~ThermodynamicOutput();

  /// writes the initial particle information of an event
  void at_eventstart(const Particles &particles,
                     const int event_number) override;

  /// writes the final particle information of an event
  void at_eventend(const Particles &particles, const int event_number,
                   double impact_parameter) override;

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
  RenamingFilePtr file_;
  // Structure that holds all the information about what to printout
  const OutputParameters out_par_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_THERMODYNAMICOUTPUT_H_
