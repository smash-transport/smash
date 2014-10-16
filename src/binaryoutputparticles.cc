/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/binaryoutputparticles.h"

#include <boost/filesystem.hpp>
#include <string>

#include <include/config.h>
#include "include/clock.h"
#include "include/configuration.h"
#include "include/forwarddeclarations.h"
#include "include/inputfunctions.h"
#include "include/particles.h"

namespace Smash {

BinaryOutputParticles::BinaryOutputParticles(bf::path path,
                                             Configuration &&config)
    : BinaryOutputBase(
          std::fopen(((path / "particles_binary.bin")).native().c_str(), "wb")),
      only_final_(config.has_value({"Only_Final"}) ? config.take({"Only_Final"})
                                                   : true) {
  /*!\Userguide
   * \page input_binary_particles Binary_particles
   * Writes the particle list at fixed times in binary format.
   *
   * \key only_final (bool, optional, default = true): \n
   * true - only final particle list at the end of each event \n
   * false - particle list output at every output interval including initial 
   * time
   * 
   * Detailed specification of the binary format can be found here:
   * \ref format_binary_
   */
  fwrite("SMSH", 4, 1, file_.get());  // magic number
  write(0);              // file format version number
  write(VERSION_MAJOR);  // version
}

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
   * 1 byte chars. Integer variables will be marked as (int), double as (d),
   * char as (char).\n
   *
   * As for OSCAR ASCII output there are two kinds of binary output:
   * particles and collisions.
   * Specifics for both particles and collisions output are the following:\n
   * **Header**
   * \code
   * magic_number, format_version(int), len(int), smash_version
   * \endcode
   * \li magic_number - 4 bytes that in ASCII read as "SMSH".
   * \li Format version is an integer number, currently it is 0.
   * \li len is the length of smash version string
   * \li smash_version is len chars that give information about SMASH version.
   *
   * **Output block header**\n
   * At start of event, end of event or any other particle output:
   * \code
   * 'p'(char) npart(int)
   * \endcode
   * \li npart is number of particle lines in the block that follows
   *
   * At interaction:
   * \code
   * 'i'(char) nin(int) nout(int)
   * \endcode
   * \li nin, nout are numbers of incoming and outgoing particles
   *
   * Block header is followed by \c nin + \c nout particle lines.
   *
   * **Particle line**
   * \code
   * t x y z mass p0 px py pz pdg(int) ID(int)
   * \endcode
   *
   * **Event end line**
   * \code
   * 'f'(char) event_number(int)
   * \endcode
   *
   * Particles output
   * ----------------
   * Written to \c particles_binary.bin file. Contains the current particle
   * list at specific moments of time. Every moment of time
   * is written as a 'p' block. For options of this output see
   * \ref input_general_, \ref input_binary_particles.
   **/

void BinaryOutputParticles::at_eventstart(const Particles &particles,
                                 const int /*event_number*/) {
  char pchar = 'p';
  if (!only_final_) {
    std::fwrite(&pchar, sizeof(char), 1, file_.get());
    write(particles.size());
    write(particles);
  }
}

void BinaryOutputParticles::at_eventend(const Particles &particles,
                               const int event_number) {
  char pchar = 'p';
  if (only_final_) {
    std::fwrite(&pchar, sizeof(char), 1, file_.get());
    write(particles.size());
    write(particles);
  }

  // Event end line
  char fchar = 'f';
  std::fwrite(&fchar, sizeof(char), 1, file_.get());
  write(event_number);

  /* Flush to disk */
  std::fflush(file_.get());
}

void BinaryOutputParticles::at_interaction(const ParticleList &/*incoming*/,
                                     const ParticleList &/*outgoing*/) {
  /* No output of this kind in particles output */
}

void BinaryOutputParticles::at_intermediate_time(const Particles &particles,
                                      const int /*event_number*/,
                                      const Clock &) {
  char pchar = 'p';
  if (!only_final_) {
    std::fwrite(&pchar, sizeof(char), 1, file_.get());
    write(particles.size());
    write(particles);
  }
}

}  // namespace Smash
