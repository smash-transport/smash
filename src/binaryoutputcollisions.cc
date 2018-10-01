/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/binaryoutputcollisions.h"

#include <string>

#include <boost/filesystem.hpp>

#include "smash/action.h"
#include "smash/clock.h"
#include "smash/config.h"
#include "smash/particles.h"

namespace smash {

BinaryOutputCollisions::BinaryOutputCollisions(const bf::path &path,
                                               std::string name,
                                               const OutputParameters &out_par)
    : BinaryOutputBase(
          path / ((name == "Collisions" ? "collisions_binary" : name) + ".bin"),
          "wb", name, out_par.get_coll_extended(name)),
      print_start_end_(out_par.coll_printstartend) {}

/*!\Userguide
 * \page format_binary_ Binary format
 *
 * Collisions output
 * -----------------
 * Written to \c collisions_binary.bin file. Contains interactions
 * (collisions, decays, box wall crossings) and optionally initial
 * and final configuration. Interactions are written in comp. frame
 * time-ordered fashion, in 'i' blocks, which includes the informations
 * of the incoming and the outgoing particles of each reaction written
 * in the 'incoming' and 'outgoing' blocks respectively.
 * Initial and final states are written as 'p' blocks. The process IDs
 * indicating the types of the reaction, such as resonance decay,
 * elastic scattering, soft string process, hard string process, etc.,
 * are written in the 'process_type' blocks. For options of this output see
 * \ref output_content_specific_options_ "content-specific output options".
 *
 * See also \ref collisions_output_in_box_modus_.
 **/

void BinaryOutputCollisions::at_eventstart(const Particles &particles,
                                           const int) {
  char pchar = 'p';
  if (print_start_end_) {
    std::fwrite(&pchar, sizeof(char), 1, file_.get());
    write(particles.size());
    write(particles);
  }
}

void BinaryOutputCollisions::at_eventend(const Particles &particles,
                                         const int event_number,
                                         double impact_parameter) {
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
  write(impact_parameter);

  // Flush to disk
  std::fflush(file_.get());
}

void BinaryOutputCollisions::at_interaction(const Action &action,
                                            const double density) {
  char ichar = 'i';
  std::fwrite(&ichar, sizeof(char), 1, file_.get());
  write(action.incoming_particles().size());
  write(action.outgoing_particles().size());
  std::fwrite(&density, sizeof(double), 1, file_.get());
  const double weight = action.get_total_weight();
  std::fwrite(&weight, sizeof(double), 1, file_.get());
  const double partial_weight = action.get_partial_weight();
  std::fwrite(&partial_weight, sizeof(double), 1, file_.get());
  const auto type = static_cast<uint32_t>(action.get_type());
  std::fwrite(&type, sizeof(uint32_t), 1, file_.get());
  write(action.incoming_particles());
  write(action.outgoing_particles());
}

BinaryOutputBase::BinaryOutputBase(const bf::path &path,
                                   const std::string &mode,
                                   const std::string &name,
                                   bool extended_format)
    : OutputInterface(name), file_{path, mode}, extended_(extended_format) {
  std::fwrite("SMSH", 4, 1, file_.get());  // magic number
  write(format_version_);                  // file format version number
  std::uint16_t format_variant = static_cast<uint16_t>(extended_);
  write(format_variant);
  write(VERSION_MAJOR);  // SMASH version
}

// write functions:
void BinaryOutputBase::write(const std::string &s) {
  const auto size = boost::numeric_cast<uint32_t>(s.size());
  std::fwrite(&size, sizeof(std::uint32_t), 1, file_.get());
  std::fwrite(s.c_str(), s.size(), 1, file_.get());
}

void BinaryOutputBase::write(const double x) {
  std::fwrite(&x, sizeof(x), 1, file_.get());
}

void BinaryOutputBase::write(const FourVector &v) {
  std::fwrite(v.begin(), sizeof(*v.begin()), 4, file_.get());
}

void BinaryOutputBase::write(const Particles &particles) {
  for (const auto &p : particles) {
    write_particledata(p);
  }
}

void BinaryOutputBase::write(const ParticleList &particles) {
  for (const auto &p : particles) {
    write_particledata(p);
  }
}

void BinaryOutputBase::write_particledata(const ParticleData &p) {
  write(p.position());
  double mass = p.effective_mass();
  std::fwrite(&mass, sizeof(mass), 1, file_.get());
  write(p.momentum());
  write(p.pdgcode().get_decimal());
  write(p.id());
  write(p.type().charge());
  if (extended_) {
    const auto history = p.get_history();
    write(history.collisions_per_particle);
    write(p.formation_time());
    write(p.xsec_scaling_factor());
    write(history.id_process);
    write(static_cast<uint32_t>(history.process_type));
    write(history.time_last_collision);
    write(history.p1.get_decimal());
    write(history.p2.get_decimal());
  }
}

}  // namespace smash
