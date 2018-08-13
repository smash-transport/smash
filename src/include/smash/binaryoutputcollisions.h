/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_BINARYOUTPUTCOLLISIONS_H_
#define SRC_INCLUDE_BINARYOUTPUTCOLLISIONS_H_

#include <string>

#include <boost/numeric/conversion/cast.hpp>

#include "file.h"
#include "forwarddeclarations.h"
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
  explicit BinaryOutputBase(const bf::path &path, const std::string &mode,
                            const std::string &name, bool extended_format);

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
  void write(const size_t x) { write(boost::numeric_cast<uint32_t>(x)); }

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
  uint16_t format_version_ = 6;
  /// Option for extended output
  bool extended_;
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
  BinaryOutputCollisions(const bf::path &path, std::string name,
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
   */
  void at_eventend(const Particles &particles, const int event_number,
                   double impact_parameter) override;

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

}  // namespace smash

#endif  // SRC_INCLUDE_BINARYOUTPUTCOLLISIONS_H_
