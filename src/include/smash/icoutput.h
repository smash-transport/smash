/*
 *
 *    Copyright (c) 2019-2020,2022-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_ICOUTPUT_H_
#define SRC_INCLUDE_SMASH_ICOUTPUT_H_

#include <filesystem>
#include <memory>
#include <string>

#include "file.h"
#include "outputinterface.h"
#include "outputparameters.h"
#include "smash/config.h"
#include "smash/outputformatter.h"
namespace smash {

/**
 * \ingroup output
 *
 * SMASH output in a format containing initial conditions for hydrodynamic
 * codes ("For_vHLLE"). Formatted such that it can be directly processed by
 * vHLLE \iref{Karpenko:2015xea}.
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
  ICOutput(const std::filesystem::path &path, const std::string &name,
           const OutputParameters &out_par);
  ~ICOutput();

  /**
   * Write event start line.
   * \param[in] event_label Numbers of the current event and ensemble.
   * \param[in] event Event info, see \ref event_info
   */
  void at_eventstart(const Particles &, const EventLabel &event_label,
                     const EventInfo &event) override;

  /**
   * Write event end line.
   * \param[in] particles Particles at end of event, expected to be empty
   * \param[in] event_label Numbers of the current event and ensemble.
   * \param[in] event Event info, see \ref event_info
   */
  void at_eventend(const Particles &particles, const EventLabel &event_label,
                   const EventInfo &event) override;

  /**
   * Unused, but needed since virtually declared in mother class.
   */
  void at_intermediate_time(const Particles &, const std::unique_ptr<Clock> &,
                            const DensityParameters &, const EventLabel &,
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
   * proper time remains unchanged during the evolution. Determined from the
   * actually removed particles.
   * By construction, tau > 0. Nevertheless it is initialized with a negative
   * number to easily find the first particle that is removed from the evolution
   * in at_interaction().
   */
  double IC_proper_time_ = -1.0;

  /// Formatter of the output
  OutputFormatter<ToASCII> formatter_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_ICOUTPUT_H_
