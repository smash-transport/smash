/*
 *
 *    Copyright (c) 2014-2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/binaryoutput.h"

#include <string>

#include <boost/filesystem.hpp>

#include "smash/action.h"
#include "smash/clock.h"
#include "smash/config.h"
#include "smash/particles.h"

namespace smash {

static constexpr int HyperSurfaceCrossing = LogArea::HyperSurfaceCrossing::id;

/*!\Userguide
 * \page format_binary_ Binary Format
 * SMASH supports a binary output version similar to the OSCAR 2013 standard.
 * It is faster to read and write and theoretically needs less disk space.
 * However, currently in ASCII OSCAR 2013 only 5 digits after the comma are
 * written for any real number, while the binary saves the whole double
 * (16 digits). By accident, this makes the sizes of the binary output files
 * approximately the same as the OSCAR ASCII files.
 * **The binary format follows the general block structure of the OSCAR
 * format:**
 * \ref oscar_general_. However, for the binary format, the data type
 * specification is stricter. The types used for the output are 4 bytes signed
 * integers, 8 bytes doubles and 1 byte chars.
 *
 * As for OSCAR ASCII output there are two kinds of binary output:
 * particles and collisions.
 * The specifics for both particles and collisions output are the following:\n
 * **Header**
 * \code
 * 4*char        uint16_t        uint16_t        uint32_t  len*char
 * magic_number, format_version, format_variant, len,      smash_version
 * \endcode
 * \li magic_number - 4 bytes that in ASCII read as "SMSH".
 * \li Format version is an integer number, currently it is 7.
 * \li Format variant is an integer number: 0 for default, 1 for extended.
 * \li len is the length of smash version string
 * \li smash_version is len chars that give information about the SMASH version.
 *
 * **Output block header**\n
 * At start of event, end of event or any other particle output:
 * \code
 * char uint32_t
 * 'p'  n_part_lines
 * \endcode
 * \li \c n_part_lines is the number of particle lines in the block that follows
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
 *        9*double          int32_t int32_t int32_t
 * t x y z mass p0 px py pz    pdg     ID      charge
 * \endcode
 *
 * **Extended Particle line**
 * <div class="fragment">
 * <div class="line">
 *   9*double       int32_t int32_t int32_t int32_t double
 *     double                    int32_t            int32_t
 *     double        int32_t         int32_t
 * </div>
 * <div class="line">
 * t x y z mass p0 px py pz pdg ID charge ncoll form_time xsecfac
 * proc_id_origin proc_type_origin time_last_coll pdg_mother1 pdg_mother2
 * </div></div>
 * \li \key proc_id_origin, \key proc_type_orgin record the id and type of
 * the last reaction that the particle has experienced.
 * \li \key time_last_coll records the time of the particle's last interaction
 * (except wall crossing), from which we can calculate the position of this
 * particle at kinetic freeze-out.
 * \li \key pdg_mother1, \key pdg_mother2 record the pdg numbers of the
 * incoming particles of the reaction where this particle is produced. If the
 * particle is produced in a resonance decay, then pdg_mother2 is set equal
 * to 0. If it is produced in a thermal bubble, then both the pdg_mother1 and
 * pdg_mother2 are set equal to zero. Both the pdg numbers are not affected
 * by elastic scatterings.
 *
 * **Event end line**
 * \code
 * char    uint32_t      double      char
 * 'f' event_number impact_parameter empty
 * \endcode
 * Where
 * \li \key event_number: Number of the event, starting with 0.
 * \li \key impact_parameter: Impact parameter [fm] of the collision in case of
 * a collider setup, 0.0 otherwise.
 * \li \key empty: 0 if there was an interaction between the projectile
 * and the target, 1 otherwise. For non-collider setups, this is always 0.
 *
 * Particles output
 * ----------------
 * The particles output is Written to the \c particles_binary.bin file.
 * It contains the current particle list at specific moments of time. Every
 * moment of time is written as a 'p' block. For options of this output see
 * \ref output_content_specific_options_ "content-specific output options".
 *
 * Collisions output
 * -----------------
 * The collisions output is Written to the \c collisions_binary.bin file.
 * It contains interactions (collisions, decays, box wall crossings) and
 * optionally the initial and final configuration. The interactions are written
 * in computational frame time-ordered fashion, in 'i' blocks, which contains
 * the information of the incoming and the outgoing particles of each reaction
 * written in the 'incoming' and 'outgoing' blocks respectively.
 * Initial and final states are written as 'p' blocks. The process IDs
 * indicating the types of the reaction, such as resonance decay,
 * elastic scattering, soft string process, hard string process, etc.,
 * are written in the 'process_type' blocks. For options of this output see
 * \ref output_content_specific_options_ "content-specific output options".
 *
 * See also \ref collisions_output_in_box_modus_.
 **/

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
void BinaryOutputBase::write(const char c) {
  std::fwrite(&c, sizeof(char), 1, file_.get());
}

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
    write(static_cast<int32_t>(history.process_type));
    write(history.time_last_collision);
    write(history.p1.get_decimal());
    write(history.p2.get_decimal());
  }
}

BinaryOutputCollisions::BinaryOutputCollisions(const bf::path &path,
                                               std::string name,
                                               const OutputParameters &out_par)
    : BinaryOutputBase(
          path / ((name == "Collisions" ? "collisions_binary" : name) + ".bin"),
          "wb", name, out_par.get_coll_extended(name)),
      print_start_end_(out_par.coll_printstartend) {}

void BinaryOutputCollisions::at_eventstart(const Particles &particles,
                                           const int, const EventInfo &) {
  const char pchar = 'p';
  if (print_start_end_) {
    std::fwrite(&pchar, sizeof(char), 1, file_.get());
    write(particles.size());
    write(particles);
  }
}

void BinaryOutputCollisions::at_eventend(const Particles &particles,
                                         const int32_t event_number,
                                         const EventInfo &event) {
  const char pchar = 'p';
  if (print_start_end_) {
    std::fwrite(&pchar, sizeof(char), 1, file_.get());
    write(particles.size());
    write(particles);
  }

  // Event end line
  const char fchar = 'f';
  std::fwrite(&fchar, sizeof(char), 1, file_.get());
  write(event_number);
  write(event.impact_parameter);
  const char empty = event.empty_event;
  write(empty);

  // Flush to disk
  std::fflush(file_.get());
}

void BinaryOutputCollisions::at_interaction(const Action &action,
                                            const double density) {
  const char ichar = 'i';
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

BinaryOutputParticles::BinaryOutputParticles(const bf::path &path,
                                             std::string name,
                                             const OutputParameters &out_par)
    : BinaryOutputBase(path / "particles_binary.bin", "wb", name,
                       out_par.part_extended),
      only_final_(out_par.part_only_final) {}

void BinaryOutputParticles::at_eventstart(const Particles &particles, const int,
                                          const EventInfo &) {
  const char pchar = 'p';
  if (only_final_ == OutputOnlyFinal::No) {
    std::fwrite(&pchar, sizeof(char), 1, file_.get());
    write(particles.size());
    write(particles);
  }
}

void BinaryOutputParticles::at_eventend(const Particles &particles,
                                        const int event_number,
                                        const EventInfo &event) {
  const char pchar = 'p';
  if (!(event.empty_event && only_final_ == OutputOnlyFinal::IfNotEmpty)) {
    std::fwrite(&pchar, sizeof(char), 1, file_.get());
    write(particles.size());
    write(particles);
  }

  // Event end line
  const char fchar = 'f';
  std::fwrite(&fchar, sizeof(char), 1, file_.get());
  write(event_number);
  write(event.impact_parameter);
  const char empty = event.empty_event;
  write(empty);

  // Flush to disk
  std::fflush(file_.get());
}

void BinaryOutputParticles::at_intermediate_time(const Particles &particles,
                                                 const std::unique_ptr<Clock> &,
                                                 const DensityParameters &,
                                                 const EventInfo &) {
  const char pchar = 'p';
  if (only_final_ == OutputOnlyFinal::No) {
    std::fwrite(&pchar, sizeof(char), 1, file_.get());
    write(particles.size());
    write(particles);
  }
}

BinaryOutputInitialConditions::BinaryOutputInitialConditions(
    const bf::path &path, std::string name, const OutputParameters &out_par)
    : BinaryOutputBase(path / "SMASH_IC.bin", "wb", name, out_par.ic_extended) {
}

void BinaryOutputInitialConditions::at_eventstart(const Particles &, const int,
                                                  const EventInfo &) {}

void BinaryOutputInitialConditions::at_eventend(const Particles &particles,
                                                const int event_number,
                                                const EventInfo &event) {
  // Event end line
  const char fchar = 'f';
  std::fwrite(&fchar, sizeof(char), 1, file_.get());
  write(event_number);
  write(event.impact_parameter);
  const char empty = event.empty_event;
  write(empty);

  // Flush to disk
  std::fflush(file_.get());

  // If the runtime is too short some particles might not yet have
  // reached the hypersurface. Warning is printed.
  if (particles.size() != 0) {
    logg[HyperSurfaceCrossing].warn(
        "End time might be too small for initial conditions output. "
        "Hypersurface has not yet been crossed by ",
        particles.size(), " particle(s).");
  }
}

void BinaryOutputInitialConditions::at_interaction(const Action &action,
                                                   const double) {
  if (action.get_type() == ProcessType::HyperSurfaceCrossing) {
    const char pchar = 'p';
    std::fwrite(&pchar, sizeof(char), 1, file_.get());
    write(action.incoming_particles().size());
    write(action.incoming_particles());
  }
}
}  // namespace smash
