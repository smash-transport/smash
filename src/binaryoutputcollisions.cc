/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/binaryoutputcollisions.h"

#include <boost/filesystem.hpp>
#include <string>

#include "include/clock.h"
#include "include/config.h"
#include "include/configuration.h"
#include "include/forwarddeclarations.h"
#include "include/inputfunctions.h"
#include "include/particles.h"

namespace Smash {

BinaryOutputCollisions::BinaryOutputCollisions(const bf::path &path,
                                               const std::string &name)
    : BinaryOutputBase(std::fopen(
          (path / (name + ".bin")).native().c_str(), "wb")),
      print_start_end_(false) {
}

BinaryOutputCollisions::BinaryOutputCollisions(const bf::path &path,
                                               Configuration &&config)
    : BinaryOutputBase(std::fopen(
          ((path / "collisions_binary.bin")).native().c_str(), "wb")),
      print_start_end_(config.take({"Print_Start_End"}, false)) {
  /*!\Userguide
   * \page input_binary_collisions Binary_collisions
   * Saves information about every collision, decay and box
   * wall crossing in a binary format. Optionally initial and
   * final particle configurations can be written out.
   *
   * \key Enable (bool, optional, default = false):\n
   * true - binary collision output enabled\n
   * false - no binary collision output
   *
   * \key Print_Start_End (bool, optional, default = false): \n
   * false - only information about collisions, decays and
   * box wall crossings during the whole evolution \n
   * true - initial and final configuration are written in addition
   *
   * Detailed specification of the binary format can be found here:
   * \ref format_binary_
   */
}

  /*!\Userguide
   * \page format_binary_ Binary format
   *
   * Collisions output
   * -----------------
   * Written to \c collisions_binary.bin file. Contains interactions
   * (collisions, decays, box wall crossings) and optionally initial
   * and final configuration. Interactions are written in comp. frame
   * time-ordered fashion, in 'i' blocks. Initial and final states
   * are written as 'p' blocks. For options of this output see
   * \ref input_general_, \ref input_binary_collisions.
   *
   * See also \ref collisions_output_in_box_modus_.
   **/


void BinaryOutputCollisions::at_eventstart(const Particles &particles,
                                 const int /*event_number*/) {
  char pchar = 'p';
  if (print_start_end_) {
    std::fwrite(&pchar, sizeof(char), 1, file_.get());
    write(particles.size());
    write(particles);
  }
}

void BinaryOutputCollisions::at_eventend(const Particles &particles,
                               const int event_number) {
  char pchar = 'p';
  if (print_start_end_) {
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

void BinaryOutputCollisions::at_interaction(const ParticleList &incoming,
                             const ParticleList &outgoing,
                             const double density,
                             const double total_cross_section,
                             const ProcessType process_type) {
  char ichar = 'i';
  std::fwrite(&ichar, sizeof(char), 1, file_.get());
  write(incoming.size());
  write(outgoing.size());
  std::fwrite(&density, sizeof(double), 1, file_.get());
  std::fwrite(&total_cross_section, sizeof(double), 1, file_.get());
  std::fwrite(&process_type, sizeof(int), 1, file_.get());
  write(incoming);
  write(outgoing);
}

void BinaryOutputCollisions::at_intermediate_time(
                                      const Particles &/*particles*/,
                                      const int /*event_number*/,
                                      const Clock &) {
  /* No output of this kind in collisions output */
}


BinaryOutputBase::BinaryOutputBase(FILE *f) : file_{f} {
  std::fwrite("SMSH", 4, 1, file_.get());  // magic number
  write(format_version_);             // file format version number
  write(VERSION_MAJOR);               // SMASH version
}


// write functions:
void BinaryOutputBase::write(const std::string &s) {
  std::int32_t size = s.size();
  std::fwrite(&size, sizeof(std::int32_t), 1, file_.get());
  std::fwrite(s.c_str(), s.size(), 1, file_.get());
}

void BinaryOutputBase::write(const FourVector &v) {
  std::fwrite(v.begin(), sizeof(*v.begin()), 4, file_.get());
}

void BinaryOutputBase::write(const Particles &particles) {
  for (const auto &p : particles) {
    //write(p.position());
    //double mass= p.pole_mass();
    //std::fwrite(&mass, sizeof(mass), 1, file_.get());
    //write(p.momentum());
    //write(p.pdgcode().get_decimal());
    //write(p.id());
    write_particledata(p);
  }
}

void BinaryOutputBase::write(const ParticleList &particles) {
  for (const auto &p : particles) {
    //write(p.momentum());
    //write(p.position());
    //write(p.pdgcode().get_decimal());
    //write(p.id());
    write_particledata(p);
  }
}

void BinaryOutputBase::write_particledata(const ParticleData &p) {
  write(p.position());
  double mass = p.pole_mass();
  std::fwrite(&mass, sizeof(mass), 1, file_.get());
  write(p.momentum());
  write(p.pdgcode().get_decimal());
  write(p.id());
} 

}  // namespace Smash
