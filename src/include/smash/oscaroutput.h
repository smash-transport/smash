/*
 *
 *    Copyright (c) 2014-2020,2022,2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_OSCAROUTPUT_H_
#define SRC_INCLUDE_SMASH_OSCAROUTPUT_H_

#include <memory>
#include <string>
#include <vector>

#include "file.h"
#include "forwarddeclarations.h"
#include "outputinterface.h"
#include "outputparameters.h"
#include "smash/outputformatter.h"

namespace smash {

/**
 * \addtogroup output
 * @{
 */

/// Selector for the output format of OscarOutput
enum OscarOutputFormat {
  OscarFormat2013,
  OscarFormat2013Extended,
  OscarFormat1999,
  ASCII
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
  OscarParticlesAtEventend = 0x008,
  /// store the state at the end of each event if it is not empty (at_eventend)
  OscarParticlesAtEventendIfNotEmpty = 0x010,
  /// store the particles that are removed on the hypersurface
  OscarParticlesIC = 0x020
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
   * \param[in] quantities List of quantities present in the output file.
   */
  OscarOutput(const std::filesystem::path &path, const std::string &name,
              const std::vector<std::string> quantities = {});

  /**
   * Writes the initial particle information of an event to the oscar output.
   * \param[in] particles Current list of all particles.
   * \param[in] event_label Numbers of event and ensemble.
   * \param[in] event Event info, see \ref event_info
   */
  void at_eventstart(const Particles &particles, const EventLabel &event_label,
                     const EventInfo &event) override;

  /**
   * Writes the final particle information of an event to the oscar output.
   * \param[in] particles Current list of particles.
   * \param[in] event_label Numbers of event and ensemble.
   * \param[in] event Event info, see \ref event_info
   */
  void at_eventend(const Particles &particles, const EventLabel &event_label,
                   const EventInfo &event) override;

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
   * \param[in] event_label Numbers of event and ensemble.
   * \param[in] event Event info, see \ref event_info
   */
  void at_intermediate_time(const Particles &particles,
                            const std::unique_ptr<Clock> &clock,
                            const DensityParameters &dens_param,
                            const EventLabel &event_label,
                            const EventInfo &event) override;

 private:
  /**
   * Write single particle information line to output.
   * \param[in] data Data of particle.
   */
  void write_particledata(const ParticleData &data);

  /**
   * Write a ToASCII::type buffer to output.
   */
  void write(const ToASCII::type &buffer);

  /**
   * Write the particle information of a list of particles to the output.
   * One line per particle.
   * \param[in] particles List of particles to be written
   */
  void write(const Particles &particles);

  /// Full filepath of the output file.
  RenamingFilePtr file_;

  /// Formatter of the output
  OutputFormatter<ToASCII> formatter_;

  /**
   * \brief Writes particle data in multiple chunks if the total buffer size
   *        exceeds a predefined maximum.
   *
   * This method avoids creating a single excessively large binary buffer when
   * writing many particles at once. Instead, it splits the write into several
   * smaller chunks. This can prevent excessive memory usage and improve
   * stability on systems or filesystems that may have trouble with very large
   * write calls.
   *
   * The maximum buffer size is currently set to 1 GB (10^9 bytes), but this can
   * be adapted in the future if needed. If the total data size of the particle
   * block is below this threshold, the method simply delegates to the regular
   * `write` function to perform a single write call.
   *
   * Otherwise, the data is accumulated particle by particle until the buffer
   * reaches the threshold. The buffer is then flushed to disk, cleared, and the
   * process continues.
   *
   * \tparam Range Container type — enforced to be either `Particles` or
   *         `ParticleList`.
   * \param[in] particles Container of particles whose particle_line
   * representation is to be written.
   */
  template <typename Converter, class Range,
            std::enable_if_t<std::is_same_v<Range, Particles> ||
                                 std::is_same_v<Range, ParticleList>,
                             bool> = true>
  void write_in_chunk(const Range &particles) {
    const auto first = std::begin(particles);
    const auto last = std::end(particles);
    if (first == last)
      return;

    constexpr std::size_t max_buffer_size = static_cast<std::size_t>(1e9);
    const std::size_t bytes_per_particle =
        formatter_.compute_single_size(*first);

    if (particles.size() <= max_buffer_size / bytes_per_particle) {
      write(formatter_.particles_chunk(particles));
      return;
    }

    typename Converter::type buffer;
    buffer.reserve(max_buffer_size);

    for (const ParticleData &particle : particles) {
      auto particle_line = formatter_.particle_line(particle);
      const std::size_t line_size = particle_line.size();
      if (buffer.size() + line_size > max_buffer_size) {
        write(buffer);
        buffer.clear();
      }
      buffer.insert(buffer.end(),
                    std::make_move_iterator(particle_line.begin()),
                    std::make_move_iterator(particle_line.end()));
    }

    if (!buffer.empty()) {
      write(buffer);
    }
  }
};

/**
 * \return A new OscarOutput object using information from \p config to
 * select the correct implementation.
 *
 * \param[in] format The output format as string, e.g. \c "Oscar2013"
 * \param[in] content The output content as string, e.g. \c "Particles"
 * \param[in] path The path to the output directory where the file(s) will be
 *             placed.
 * \param[in] out_par A structure containing parameters of the output, in
 *             particular if it is extended or not, if printing only final
 *             particles in event, etc.
 */
std::unique_ptr<OutputInterface> create_oscar_output(
    const std::string &format, const std::string &content,
    const std::filesystem::path &path, const OutputParameters &out_par);

// @}

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_OSCAROUTPUT_H_
