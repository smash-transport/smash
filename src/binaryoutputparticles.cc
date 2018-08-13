/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/binaryoutputparticles.h"

#include <string>

#include <boost/filesystem.hpp>

#include "smash/clock.h"
#include "smash/config.h"
#include "smash/forwarddeclarations.h"
#include "smash/particles.h"

namespace smash {

BinaryOutputParticles::BinaryOutputParticles(const bf::path &path,
                                             std::string name,
                                             const OutputParameters &out_par)
    : BinaryOutputBase(path / "particles_binary.bin", "wb", name,
                       out_par.part_extended),
      only_final_(out_par.part_only_final) {}

/*!\Userguide
 * \page format_binary_ Binary format
 * SMASH supports a binary version of output similar to OSCAR 2013 standard.
 * It is faster to read and write and theoretically needs less disk space.
 * However, currently in ASCII OSCAR 2013 only 5 digits after comma are
 * written for any real number, while binary saves the whole double
 * (16 digits). By accident this makes sizes of binary output files
 * approximately the same as OSCAR ASCII files.
 * **The format follows general block structure of OSCAR format:**
 * \ref oscar_general_. However, for binary specification is stricter.
 * Types used for output are 4 bytes signed integers, 8 bytes doubles and
 * 1 byte chars.
 *
 * As for OSCAR ASCII output there are two kinds of binary output:
 * particles and collisions.
 * Specifics for both particles and collisions output are the following:\n
 * **Header**
 * \code
 * 4*char        uint16_t        uint16_t        uint32_t  len*char
 * magic_number, format_version, format_variant, len,      smash_version
 * \endcode
 * \li magic_number - 4 bytes that in ASCII read as "SMSH".
 * \li Format version is an integer number, currently it is 4.
 * \li Format variant is an integer number: 0 for default, 1 for extended.
 * \li len is the length of smash version string
 * \li smash_version is len chars that give information about SMASH version.
 *
 * **Output block header**\n
 * At start of event, end of event or any other particle output:
 * \code
 * char uint32_t
 * 'p'  npart
 * \endcode
 * \li \c npart is number of particle lines in the block that follows
 *
 * At interaction:
 * \code
 * char uint32_t uint32_t double  double   uint32_t
 * 'i'  nin      nout     density xsection process_type
 * \endcode
 * \li \c nin, \c nout are numbers of incoming and outgoing particles
 *
 * Block header is followed by \c nin + \c nout particle lines.
 *
 * **Particle line**
 * \code
 *     9*double             int int int
 * t x y z mass p0 px py pz pdg ID charge
 * \endcode
 *
 * **Extended Particle line**
 * <div class="fragment">
 * <div class="line">
 *   9*double        int int int int     double
 *     double                    uint32_t            uint32_t
 *     double        int         int
 * </div>
 * <div class="line">
 * t x y z mass p0 px py pz pdg ID charge Ncoll formation_time
 * xs_scaling_factor process_ID_origin process_type_origin
 * time_of_origin PDG_mother1 PDG_mother2
 * </div></div>
 *
 * **Event end line**
 * \code
 * char    int            int
 * 'f' event_number impact_parameter
 * \endcode
 *
 * Particles output
 * ----------------
 * Written to \c particles_binary.bin file. Contains the current particle
 * list at specific moments of time. Every moment of time
 * is written as a 'p' block. For options of this output see
 * \ref output_content_specific_options_ "content-specific output options".
 **/

void BinaryOutputParticles::at_eventstart(const Particles &particles,
                                          const int) {
  char pchar = 'p';
  if (!only_final_) {
    std::fwrite(&pchar, sizeof(char), 1, file_.get());
    write(particles.size());
    write(particles);
  }
}

void BinaryOutputParticles::at_eventend(const Particles &particles,
                                        const int event_number,
                                        double impact_parameter) {
  char pchar = 'p';
  /// \todo(schaefer) Is the if-statement below really correct, if only_final_
  /// is False?
  if (only_final_) {
    std::fwrite(&pchar, sizeof(char), 1, file_.get());
    write(particles.size());
    write(particles);
  }

  // Event end line
  char fchar = 'f';
  std::fwrite(&fchar, sizeof(char), 1, file_.get());
  write(event_number);
  write(impact_parameter);

  // Flush to disk
  std::fflush(file_.get());
}

void BinaryOutputParticles::at_intermediate_time(const Particles &particles,
                                                 const Clock &,
                                                 const DensityParameters &) {
  char pchar = 'p';
  if (!only_final_) {
    std::fwrite(&pchar, sizeof(char), 1, file_.get());
    write(particles.size());
    write(particles);
  }
}

}  // namespace smash
