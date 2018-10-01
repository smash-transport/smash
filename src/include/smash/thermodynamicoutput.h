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
  /**
   * Construct Output
   * param[in] path Path to output
   * param[in] name Filename
   * param[in] out_par Parameters of output
   */
  ThermodynamicOutput(const bf::path &path, const std::string &name,
                      const OutputParameters &out_par);
  /// Default destructor
  ~ThermodynamicOutput();

  /**
   *  writes the event header
   *  \param[in] particles Dummy, is just here to satisfy inheritance
   *  \param[in] event_number Current event number,
   *             that will be written to the header
   */
  void at_eventstart(const Particles &particles,
                     const int event_number) override;

  /**
   *  only flushes the output the file
   *  \param[in] particles Dummy, is just here to satisfy inheritance
   *  \param[in] particles Dummy, is just here to satisfy inheritance
   *  \param[in] particles Dummy, is just here to satisfy inheritance
   */
  void at_eventend(const Particles &particles, const int event_number,
                   double impact_parameter) override;

  /**
   *  Writes thermodynamics every fixed time interval. For configuring
   *  the output see \ref thermodyn_output_user_guide_.
   *  \param[in] particles, from which thermodynamic variables are computed
   *  \param[in] clock Clock, needed to get current time
   *  \param[in] dens_param set of parameters, defining smearing. For more
   *             info about smearing see \ref thermodyn_output_user_guide_.
   */
  void at_intermediate_time(const Particles &particles, const Clock &clock,
                            const DensityParameters &dens_param) override;

  /**
   * Prints density along the specified line. Useful to make 1D plots of
   * density profiles.
   * \param[in] file_name name of the file to print out
   * \param[in] plist particles, from which density is computed
   * \param[in] dens_type type of density, see ::DensityType
   * \param[in] line_start starting point of the line
   * \param[in] line_end ending point of the line
   * \param[in] n_points number of points along the line, where density
   *            is printed out
   */
  void density_along_line(const char *file_name, const ParticleList &plist,
                          const DensityParameters &param, DensityType dens_type,
                          const ThreeVector &line_start,
                          const ThreeVector &line_end, int n_points);

 private:
  /// Pointer to output file
  RenamingFilePtr file_;
  /// Structure that holds all the information about what to printout
  const OutputParameters out_par_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_THERMODYNAMICOUTPUT_H_
