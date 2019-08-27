/*
 *
 *    Copyright (c) 2014-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_ICOUTPUT_H_
#define SRC_INCLUDE_ICOUTPUT_H_

#include <string>

#include <boost/filesystem.hpp>

#include "file.h"
#include "outputinterface.h"
#include "outputparameters.h"
#include "smash/config.h"

namespace smash {

/**
 * \ingroup output
 * SMASH output in ASCII format containing initial conditions for hydrodynamic
 * codes.
 */
class ICOutput : public OutputInterface {
 public:
  /**
   * Create a new IC output.
   *
   * \param[in] path Path to the output file.
   * \param[in] name Name of the output.
   * \param[in] out_par Additional information on the configured output.
   */
  ICOutput(const bf::path &path, const std::string &name,
           const OutputParameters &out_par);
  ~ICOutput();

  /**
   * Write event start line.
   * \param[in] particles Unused, needed since inherited.
   * \param[in] event_number Number of the current event.
   */

  void at_eventstart(const Particles &particles,
                     const int event_number) override;

  /**
   * Write event end line.
   * \param[in] particles Unused, needed since inherited.
   * \param[in] event_number Number of the current event.
   * \param[in] impact_parameter Unused, needed since inherited.
   * \param[in] empty_event Unused, needed since inherited.
   */
  void at_eventend(const Particles &particles, const int event_number,
                   double impact_parameter, bool empty_event) override;

  /**
   * Unused, but needed since virtually declared in mother class.
   * \param[in] particles Unused, needed since inherited.
   * \param[in] clock Unused, needed since inherited.
   * \param[in] dens_param Unused, needed since inherited.
   */
  void at_intermediate_time(const Particles &particles, const Clock &clock,
                            const DensityParameters &dens_param) override;
  /**
   * Write particle data at the hypersurface crossing point to the IC output.
   *
   * \param[in] particles The particle that is removed.
   * \param[in] clock Time at which hypersurface is crossed.
   * \param[in] dens_param Unused, needed since inherited.
   */
  void at_interaction(const Action &action, const double density) override;

 private:
  /// Pointer to output file
  RenamingFilePtr file_;
  /// Structure that holds all the information about what to printout
  const OutputParameters out_par_;

  /**
   * Proper time at which hypersurface is created.
   * By construction, tau > 0. Nevertheless it is initialized with a negative
   * number to easily find the first particle that is removed from the evolution
   * in at_interaction().
   */
  double IC_proper_time_ = -1.0;
};

}  // namespace smash

#endif  // SRC_INCLUDE_ICOUTPUT_H_
