
/*
 *
 *    Copyright (c) 2014-2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_HEPMCOUTPUT_H_
#define SRC_INCLUDE_SMASH_HEPMCOUTPUT_H_

#include <HepMC3/WriterAscii.h>

#include <boost/filesystem.hpp>

#include "hepmcinterface.h"

namespace smash {

/**
 * \ingroup output
 * \brief SMASH output to HepMC file
 *
 * This class writes a vertex connecting all intial particles with all final
 * particles into a HepMC outputfile. In collider mode, projectile and target
 * are combined into single intial particles with a nuclear pdg code. The output
 * file is a human-readable ASCII file. HepMC version 3 is used.
 *
 * More details of the output format can be found in the User Guide.
 */
class HepMcOutput : public HepMcInterface {
 public:
  using Base = HepMcInterface;
  /**
   * Create HepMC particle output.
   *
   * \param[in] path Output path.
   * \param[in] name Name of the output.
   * \param[in] out_par Unused, needed since inhertied.
   * \param[in] total_N Total number of particles in both nuclei.
   * \param[in] proj_N  Number of particles in projectile.
   */
  HepMcOutput(const bf::path &path, std::string name,
              const OutputParameters &out_par, const int total_N,
              const int proj_N);

  /// Destructor renames file
  ~HepMcOutput();
  /**
   * Add the final particles information of an event to the central vertex.
   * Store impact paramter and write event.
   *
   * \param[in] particles Current list of particles.
   * \param[in] event_number Number of event.
   * \param[in] event Event info, see \ref event_info
   */
  void at_eventend(const Particles &particles, const int32_t event_number,
                   const EventInfo &event) override;

 private:
  /// Filename of output
  const bf::path filename_;
  /// Filename of output as long as simulation is still running.
  bf::path filename_unfinished_;
  /// Pointer to Ascii HepMC3 output file
  std::unique_ptr<HepMC3::WriterAscii> output_file_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_HEPMCOUTPUT_H_
