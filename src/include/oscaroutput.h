/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_OSCAROUTPUT_H_
#define SRC_INCLUDE_OSCAROUTPUT_H_

#include <memory>
#include <string>

#include "filedeleter.h"
#include "forwarddeclarations.h"
#include "outputinterface.h"
#include "outputparameters.h"

namespace Smash {

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
  OscarOutput(const bf::path &path, std::string name);

  /// writes the initial particle information of an event
  void at_eventstart(const Particles &particles,
                     const int event_number) override;

  /// writes the final particle information of an event
  void at_eventend(const Particles &particles, const int event_number,
                   double impact_parameter) override;

  /// Write a prefix line and a line per particle to OSCAR output.
  void at_interaction(const Action &action, const double density) override;

  void at_intermediate_time(const Particles &particle, const Clock &clock,
                            const DensityParameters &dens_param) override;

 private:
  void write_particledata(const ParticleData &data);
  void write(const Particles &particles);

  int current_event_ = 0;
  FilePtr file_;
};

/**
 * Returns a new OscarOutput object using information from \p config to
 * select the correct implementation.
 *
 * \param format A string: "Oscar2013" or "Oscar1999"
 * \param content A string: "Particles", "Collisions", "Photons"
               or "Dileptons".
 * \param path The path to the output directory where the file(s) will be
 *             placed.
 * \param out_par A structure containing parameters of the output, in particular
 *             if it is extended or not, if printing only final particles
 *             in event, etc.
 */
std::unique_ptr<OutputInterface> create_oscar_output(
    std::string format, std::string content, const bf::path &path,
    const OutputParameters &out_par);

// @}

}  // namespace Smash

#endif  // SRC_INCLUDE_OSCAROUTPUT_H_
