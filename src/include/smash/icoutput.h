/*
 *
 *    Copyright (c) 2019-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_ICOUTPUT_H_
#define SRC_INCLUDE_SMASH_ICOUTPUT_H_

#include <string>

#include <boost/filesystem.hpp>

#include "file.h"
#include "outputinterface.h"
#include "outputparameters.h"
#include "smash/config.h"

namespace smash {

/**
 * \ingroup output
 *
 * SMASH output in ASCII format containing initial conditions for hydrodynamic
 * codes. Formatted such that it can be directly processed by vHLLE
 * \iref{Karpenko:2015xea}.
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
   * \param[in] event_number Number of the current event.
   */
  void at_eventstart(const Particles &, const int event_number,
                     const EventInfo &) override;

  /**
   * Write event end line.
   * \param[in] particles Particles at end of event, expected to be empty
   * \param[in] event_number Number of the current event.
   */
  void at_eventend(const Particles &particles, const int event_number,
                   const EventInfo &) override;

  /**
   * Unused, but needed since virtually declared in mother class.
   */
  void at_intermediate_time(const Particles &, const std::unique_ptr<Clock> &,
                            const DensityParameters &,
                            const EventInfo &) override;
  /**
   * Write particle data at the hypersurface crossing point to the IC output.
   *
   * \param[in] action Details about the action
   */
  void at_interaction(const Action &action, const double) override;

 private:
  /// Pointer to output file
  RenamingFilePtr file_;
  /// Structure that holds all the information about what to printout
  const OutputParameters out_par_;

  /**
   * Proper time of the particles removed when extracting initial conditions.
   * Parameter used for testing purposes only. Used to verify that the initial
   * proper time remains unchanged during the evolution. Dewtermined from the
   * actually removed particles.
   * By construction, tau > 0. Nevertheless it is initialized with a negative
   * number to easily find the first particle that is removed from the evolution
   * in at_interaction().
   */
  double IC_proper_time_ = -1.0;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_ICOUTPUT_H_
