/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_BINARYOUTPUTCOLLISIONS_H_
#define SRC_INCLUDE_BINARYOUTPUTCOLLISIONS_H_

#include <string>

#include "outputinterface.h"
#include "filedeleter.h"
#include "forwarddeclarations.h"
#include "configuration.h"

namespace Smash {

/**
 * \ingroup output
 * Base class for SMASH binary output.
 **/
class BinaryOutputBase : public OutputInterface {
 protected:
  BinaryOutputBase(FILE *f);
  void write(const std::string &s);
  void write(const FourVector &v);
  void write(std::int32_t x) {
    std::fwrite(&x, sizeof(x), 1, file_.get());
  }
  void write(const Particles &particles);
  void write(const ParticleList &particles);

  /// Binary particles output
  FilePtr file_;

 private:
  /// file format version number
  const int format_version_ = 2; 
};

/**
 * \ingroup output
 * \brief Saves SMASH collision history to binary file.
 *
 * This class writes every collision, decay and box wall crossing
 * to the output file. Optionally one can also write
 * initial and final particle lists to the same file.
 * Output file is binary and has a block structure.
 *
 * Details of the output format can be found
 * on the wiki in User Guide section, look for binary output.
 */
class BinaryOutputCollisions : public BinaryOutputBase {
 public:
  BinaryOutputCollisions(bf::path path, Configuration&& config);

  /// writes the initial particle information of an event
  void at_eventstart(const Particles &particles,
                     const int event_number) override;

  /// writes the final particle information of an event
  void at_eventend(const Particles &particles, const int event_number) override;

  void at_interaction(const ParticleList &incoming_particles,
                      const ParticleList &outgoing_particles,
                      const double density,
                      const double total_cross_section) override;
  /// writes particles every time interval fixed by option OUTPUT_INTERVAL
  void at_intermediate_time(const Particles &particles, const int event_number,
                          const Clock &clock) override;

 private:
  /// Option: print initial and final particles or not
  bool print_start_end_;
};


}  // namespace Smash

#endif  // SRC_INCLUDE_BINARYOUTPUTCOLLISIONS_H_
