/*
 *
 *    Copyright (c) 2014-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_HEPMCOUTPUT_H_
#define SRC_INCLUDE_HEPMCOUTPUT_H_


#include "outputinterface.h"
#include "outputparameters.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenVertex.h"


namespace smash {


const int status_code_for_beam_particles = 4;
const int status_code_for_final_particles = 1;

class HepMcOutput : public OutputInterface {
 public:
  /**
   * Create binary particle output.
   *
   * \param[in] path Output path.
   * \param[in] name Name of the output.
   * \param[in] out_par A structure containing parameters of the output.
   */
  HepMcOutput(const bf::path &path, std::string name,
              const OutputParameters &out_par);

  /**
   * Writes the initial particle information list of an event to the binary
   * output.
   * \param[in] particles Current list of all particles.
   * \param[in] event_number Unused, needed since inherited.
   */
  void at_eventstart(const Particles &particles,
                     const int event_number) override;

  /**
   * Writes the final particle information list of an event to the binary
   * output.
   * \param[in] particles Current list of particles.
   * \param[in] event_number Number of event.
   * \param[in] impact_parameter Impact parameter of this event.
   * \param[in] empty_event Whether there was no collision between target
   *            and projectile
   */
  void at_eventend(const Particles &particles, const int32_t event_number,
                   double impact_parameter, bool empty_event) override;

  /**
   * Writes an interaction block, including information about the incoming and
   * outgoing particles, to the binary output.
   * \param[in] action Action that holds the information of the interaction.
   * \param[in] density Density at the interaction point.
   */
  void at_interaction(const Action &action, const double density) override;

private:

  HepMC3::WriterAscii output_file_;
  //
  HepMC3::GenEvent current_event_;
  //
  HepMC3::GenVertexPtr vertex_;

};

}  // namespace smash

#endif  // SRC_INCLUDE_HEPMCOUTPUT_H_
