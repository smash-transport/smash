/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_OSCAROUTPUT_H_
#define SRC_INCLUDE_OSCAROUTPUT_H_

#include <memory>
#include <string>

#include "file.h"
#include "forwarddeclarations.h"
#include "outputinterface.h"
#include "outputparameters.h"

namespace smash {

/**
 * \addtogroup output
 * @{
 */

/// Selector for the output format of OscarOutput
enum OscarOutputFormat {
  OscarFormat2013,
  OscarFormat2013Extended,
  OscarFormat1999
};

/**
 * \brief Flags for the \p Contents template parameter of OscarOutput.
 *
 * Flags can be combined with binary OR operators to some arbitrary int. That's
 * why the values of the enumerators are written out (in hexadecimal), to ensure
 * every flag occupies a single bit.
 */
enum OscarOutputContents {
  /// store interaction information (write_interaction)
  OscarInteractions = 0x001,
  /// store the state after N timesteps (after_Nth_timestep)
  OscarTimesteps = 0x002,
  /// store the state at the start of each event (at_eventstart)
  OscarAtEventstart = 0x004,
  /// store the state at the end of each event (at_eventend)
  OscarParticlesAtEventend = 0x008
};

/**
 * \tparam Format Determines the variant of OSCAR formatting that is used. See
 *                OscarOutputFormat.
 * \tparam Contents Determines what information will be written to file. This
 *                  integer is a bitflag that can be constructed from ORing
 *                  enumerators from OscarOutputContents together.
 */
template <OscarOutputFormat Format, int Contents>
class OscarOutput : public OutputInterface {
 public:
  /**
   * Create oscar output.
   *
   * \param[in] path Output path.
   * \param[in] name Name of the ouput.
   */
  OscarOutput(const bf::path &path, const std::string &name);

  /**
   * Writes the initial particle information of an event to the oscar output.
   * \param[in] particles Current list of all particles.
   * \param[in] event_number Number of event.
   */
  void at_eventstart(const Particles &particles,
                     const int event_number) override;

  /**
   * Writes the final particle information of an event to the oscar output.
   * \param[in] particles Current list of particles.
   * \param[in] event_number Number of event.
   * \param[in] impact_parameter Impact parameter of this event.
   */
  void at_eventend(const Particles &particles, const int event_number,
                   double impact_parameter) override;

  /**
   * Writes a interaction prefix line and a line for every incoming and
   * outgoing particle to the oscar output.
   * \param[in] action Action that holds the information of the interaction.
   * \param[in] density Density at the interaction point.
   */
  void at_interaction(const Action &action, const double density) override;

  /**
   * Writes a prefix line then write out all current particles.
   *
   * \param[in] particles Current list of particles.
   * \param[in] clock Unused, needed since inherited.
   * \param[in] dens_param Unused, needed since inherited.
   */
  void at_intermediate_time(const Particles &particles, const Clock &clock,
                            const DensityParameters &dens_param) override;

 private:
  /**
   * Write single particle information line to output.
   * \param[in] data Data of particle.
   */
  void write_particledata(const ParticleData &data);

  /**
   * Write the particle information of a list of particles to the output.
   * One line per particle.
   * \param[in] particles List of particles to be written
   */
  void write(const Particles &particles);

  /// Keep track of event number.
  int current_event_ = 0;

  /// Full filepath of the output file.
  RenamingFilePtr file_;
};

/**
 * \return A new OscarOutput object using information from \p config to
 * select the correct implementation.
 *
 * \param[in] format A string: "Oscar2013" or "Oscar1999"
 * \param[in] content A string: "Particles", "Collisions", "Photons"
               or "Dileptons".
 * \param[in] path The path to the output directory where the file(s) will be
 *             placed.
 * \param[in] out_par A structure containing parameters of the output, in
 *             particular if it is extended or not, if printing only final 
 *             particles in event, etc.
 */
std::unique_ptr<OutputInterface> create_oscar_output(
    const std::string &format, const std::string &content, const bf::path &path,
    const OutputParameters &out_par);

// @}

}  // namespace smash

#endif  // SRC_INCLUDE_OSCAROUTPUT_H_
