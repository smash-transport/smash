/*
 *
 *    Copyright (c) 2014-2020,2022-2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_BINARYOUTPUT_H_
#define SRC_INCLUDE_SMASH_BINARYOUTPUT_H_
#include <cstdint>
#include <memory>
#include <string>

#include "file.h"
#include "forwarddeclarations.h"
#include "numeric_cast.h"
#include "outputformatter.h"
#include "outputinterface.h"
#include "outputparameters.h"

namespace smash {

/**
 * \ingroup output
 * Base class for SMASH binary output.
 */
class BinaryOutputBase : public OutputInterface {
 protected:
  /**
   * Create binary output base.
   *
   * \param[in] path Output path.
   * \param[in] mode Is used to determine the file access mode.
   * \param[in] name Name of the output.
   * \param[in] extended_format Is the written output extended.
   */
  explicit BinaryOutputBase(const std::filesystem::path &path,
                            const std::string &mode, const std::string &name,
                            bool extended_format,
                            const std::vector<std::string> &quantities = {}

  );

  /**
   * Write several bytes to the binary output.
   * \param[in] chunk vector of bytes to be written.
   *
   */
  void write(const std::vector<char> &chunk);

  /**
   * Write byte to binary output.
   * \param[in] c Value to be written.
   */
  void write(const char c);

  /**
   * Write string to binary output.
   * \param[in] s String to be written.
   */
  void write(const std::string &s);

  /**
   * Write double to binary output.
   * \param[in] x Value to be written.
   */
  void write(const double x);

  /**
   * Write four-vector to binary output.
   * \param[in] v Four-vector to be written.
   */
  void write(const FourVector &v);

  /**
   * Write integer (32 bit) to binary output.
   * \param[in] x Value to be written.
   */
  void write(const std::int32_t x) {
    std::fwrite(&x, sizeof(x), 1, file_.get());
  }

  /**
   * Write unsigned integer (32 bit) to binary output.
   * \param[in] x Value to be written.
   */
  void write(const std::uint32_t x) {
    std::fwrite(&x, sizeof(x), 1, file_.get());
  }

  /**
   * Write unsigned integer (16 bit) to binary output.
   * \param[in] x Value to be written.
   */
  void write(const std::uint16_t x) {
    std::fwrite(&x, sizeof(x), 1, file_.get());
  }

  /**
   * Write a std::size_t to binary output.
   * \param[in] x Value to be written.
   */
  void write(const size_t x) { write(smash::numeric_cast<uint32_t>(x)); }

  /**
   * Write particle data of each particle in particles to binary output.
   * \param[in] particles List of particles, whose data is to be written.
   */
  void write(const Particles &particles);

  /**
   * Write each particle data entry to binary output.
   * \param[in] particles List of particles, whose data is to be written.
   */
  void write(const ParticleList &particles);

  /**
   * Write particle data to binary output.
   * \param[in] p Particle data to be written.
   */
  void write_particledata(const ParticleData &p);

  /// Binary particles output file path
  RenamingFilePtr file_;

 private:
  /// Binary file format version number
  const uint16_t format_version_ = 10;
  /// Option for extended output
  bool extended_;
  /// If custom is used extended field is set to:
  const uint16_t format_custom = 2;

  OutputFormatter<ToBinary> formatter_;
};

/**
 * \ingroup output
 * \brief Saves SMASH collision history to binary file.
 *
 * This class writes each collision, decay and box wall crossing
 * to the output file. Optionally, one can also write the
 * initial and final particle lists to the same file.
 * The output file is binary and has a block structure.
 *
 * Details of the output format can be found
 * on the wiki in the User Guide section, look for binary output.
 */
class BinaryOutputCollisions : public BinaryOutputBase {
 public:
  /**
   * Create binary particle output.
   *
   * \param[in] path Output path.
   * \param[in] name Name of the output.
   * \param[in] out_par A structure containing parameters of the output.
   */
  BinaryOutputCollisions(const std::filesystem::path &path, std::string name,
                         const OutputParameters &out_par);

  /**
   * Writes the initial particle information list of an event to the binary
   * output.
   * \param[in] particles Current list of all particles.
   * \param[in] event_number Unused, needed since inherited.
   * \param[in] event Event info, see \ref event_info
   */
  void at_eventstart(const Particles &particles, const int event_number,
                     const EventInfo &event) override;

  /**
   * Writes the final particle information list of an event to the binary
   * output.
   * \param[in] particles Current list of particles.
   * \param[in] event_number Number of event.
   * \param[in] event Event info, see \ref event_info
   */
  void at_eventend(const Particles &particles, const int32_t event_number,
                   const EventInfo &event) override;

  /**
   * Writes an interaction block, including information about the incoming and
   * outgoing particles, to the binary output.
   * \param[in] action Action that holds the information of the interaction.
   * \param[in] density Density at the interaction point.
   */
  void at_interaction(const Action &action, const double density) override;

 private:
  /// Write initial and final particles additonally to collisions?
  bool print_start_end_;
};

/**
 * \ingroup output
 *
 * \brief Writes the particle list at specific times to the binary file
 *
 * This class writes the current particle list at a specific time t
 * to the binary output file. This specific time can
 * be: event start, event end or every next time interval \f$\Delta t \f$.
 * Writing (or not writing) the output at these moments is controlled by
 * different options. The time interval \f$\Delta t \f$ is also regulated by an
 * option. The output file is binary and has a block structure.
 *
 * Details of the output format can be found
 * on the wiki in the User Guide section, look for binary output.
 */
class BinaryOutputParticles : public BinaryOutputBase {
 public:
  /**
   * Create binary particle output.
   *
   * \param[in] path Output path.
   * \param[in] name Name of the ouput.
   * \param[in] out_par A structure containing the parameters of the output.
   */
  BinaryOutputParticles(const std::filesystem::path &path, std::string name,
                        const OutputParameters &out_par);

  /**
   * Writes the initial particle information of an event to the binary output.
   * \param[in] particles Current list of all particles.
   * \param[in] event_number Unused, needed since inherited.
   * \param[in] event Event info, see \ref event_info
   */
  void at_eventstart(const Particles &particles, const int event_number,
                     const EventInfo &event) override;

  /**
   * Writes the final particle information of an event to the binary output.
   * \param[in] particles Current list of particles.
   * \param[in] event_number Number of event.
   * \param[in] event Event info, see \ref event_info
   */
  void at_eventend(const Particles &particles, const int event_number,
                   const EventInfo &event) override;

  /**
   * Writes particles at each time interval; fixed by option OUTPUT_INTERVAL.
   * \param[in] particles Current list of particles.
   * \param[in] clock Unused, needed since inherited.
   * \param[in] dens_param Unused, needed since inherited.
   * \param[in] event Event info, see \ref event_info.
   */
  void at_intermediate_time(const Particles &particles,
                            const std::unique_ptr<Clock> &clock,
                            const DensityParameters &dens_param,
                            const EventInfo &event) override;

 private:
  /// Whether final- or initial-state particles should be written.
  OutputOnlyFinal only_final_;
};

/**
 * \ingroup output
 *
 * \brief Writes the particles when crossing the hypersurface to the binary file
 *
 * This class writes each particle to the binary output at the time of crossing
 * the hypersurface. This time corresponds to the proper time of the
 * hypersruface, which is - if not specified differently in the configuration -
 * the passing time of the two nuclei.
 *
 * Details of the output format can be found
 * on the wiki in the User Guide section, look for Output: Initial Conditions.
 */
class BinaryOutputInitialConditions : public BinaryOutputBase {
 public:
  /**
   * Create binary initial conditions particle output.
   *
   * \param[in] path Output path.
   * \param[in] name Name of the ouput.
   * \param[in] out_par A structure containing the parameters of the output.
   */
  BinaryOutputInitialConditions(const std::filesystem::path &path,
                                std::string name,
                                const OutputParameters &out_par);

  /**
   * Writes the initial particle information of an event to the binary output.
   * Function unused for IC output. Needed since inherited.
   */
  void at_eventstart(const Particles &, const int, const EventInfo &) override;

  /**
   * Writes the final particle information of an event to the binary output.
   * \param[in] particles Current list of particles.
   * \param[in] event_number Number of event.
   * \param[in] event Event info, see \ref event_info
   */
  void at_eventend(const Particles &particles, const int event_number,
                   const EventInfo &event) override;

  /**
   * Writes particles that are removed when crossing the hypersurface to the
   * output. Note that the particle information is written as a particle block,
   * not as an interaction block.
   * \param[in] action Action that holds the information of the interaction.
   */
  void at_interaction(const Action &action, const double) override;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_BINARYOUTPUT_H_
