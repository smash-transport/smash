
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

#include <memory>
#include <string>

#include "forwarddeclarations.h"
#include "outputinterface.h"
#include "outputparameters.h"

#include <boost/filesystem.hpp>

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/WriterAscii.h"

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
class HepMcOutput : public OutputInterface {
 public:
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
   * Add the initial particles information of an event to the central vertex.
   * Construct projectile and target particles with nuclear pdg code if
   * collider.
   * \param[in] particles Current list of all particles.
   * \param[in] event_number Current event number
   * \throw std::invalid_argument if nuclei with non-nucleon particle
   *                              (like hypernuclei) are tried to be constructed
   */
  void at_eventstart(const Particles &particles, const int event_number,
                     const EventInfo &) override;

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
  /// HepMC convention: status code for target and projecticle particles
  static const int status_code_for_beam_particles;
  /// HepMC convention: status code for final particles
  static const int status_code_for_final_particles;

  /// Filename of output
  const bf::path filename_;
  /// Filename of output as long as simulation is still running.
  bf::path filename_unfinished_;

  /**
   * Total number of nucleons in projectile and target,
   * needed for converting nuclei to single particles.
   */
  const int total_N_;
  /**
   * Total number of nucleons in projectile,
   * needed for converting nuclei to single particles.
   */
  const int proj_N_;

  /**
   * Construct nulcear pdg code for porjectile and target, see PDG chapter
   * "Monte Carlo particle numbering scheme" for details.
   *
   * \param[in] na Number of all particles in nucleus (A)
   * \param[in] nz Number of all charged particles in nucleus (Z).
   *
   * Note: Isomers and hypernuclei are not supported. Also note that one would
   * usually create a PdgCode object to store the PDG code of the nuclei.
   * This PdgCode object however would be first constructed with the integer
   * determined here and then directly returned again with a method of the
   * PdgCode class, making the creation of a PdgCode class instance overkill.
   * A possible improvement would be the inclusion of this function in the
   * PdgCode class.
   */
  int construct_nuclear_pdg_code(int na, int nz) const;

  /// Pointer to Ascii HepMC3 output file
  std::unique_ptr<HepMC3::WriterAscii> output_file_;

  /// HepMC3::GenEvent pointer for current event
  std::unique_ptr<HepMC3::GenEvent> current_event_;

  /// HepMC3::GenVertex pointer for central vertex in event
  HepMC3::GenVertexPtr vertex_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_HEPMCOUTPUT_H_
