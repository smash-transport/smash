
/*
 *
 *    Copyright (c) 2020-2022,2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_HEPMCOUTPUT_H_
#define SRC_INCLUDE_SMASH_HEPMCOUTPUT_H_

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "HepMC3/Writer.h"

#include "hepmcinterface.h"

namespace smash {

/**
 * \ingroup output
 * \brief SMASH output to HepMC file
 *
 * This class writes a vertex connecting all intial particles with all final
 * particles into a HepMC outputfile. In collider mode, projectile and target
 * are combined into single intial particles with a nuclear pdg code. The output
 * file can be a human-readable ASCII file or a ROOT Tree binary file.
 * HepMC version 3 is used.
 *
 * More details of the output format can be found in the User Guide.
 */
class HepMcOutput : public HepMcInterface {
 public:
  /**
   * Create HepMC particle output.
   *
   * \param[in] path Output path.
   * \param[in] name Name of the output.
   * \param[in] full_event Whether the full event or only final-state particles
                           are printed in the output
   * \param[in] HepMC3_output_type: "root" or "asciiv3"
   */
  HepMcOutput(const std::filesystem::path &path, std::string name,
              const bool full_event, std::string HepMC3_output_type);

  /// Destructor renames file
  ~HepMcOutput();
  /**
   * Add the final particles information of an event to the central vertex.
   * Store impact parameter and write event.
   *
   * \param[in] particles Current list of particles.
   * \param[in] event_label Event/ensemble numbers
   * \param[in] event Event info, see \ref event_info
   */
  void at_eventend(const Particles &particles, const EventLabel &event_label,
                   const EventInfo &event) override;

 private:
  /// Filename of output
  const std::filesystem::path filename_;
  /// Filename of output as long as simulation is still running.
  std::filesystem::path filename_unfinished_;
  /// Pointers to the base class of HepMC3 output files
  std::unique_ptr<HepMC3::Writer> output_file_;
  /// enum to identify the HepMC3 output type
  typedef enum enum_output { asciiv3, treeroot } type_of_HepMC3_output;
  /// HepMC3 output type
  type_of_HepMC3_output output_type_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_HEPMCOUTPUT_H_
