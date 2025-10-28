/*
 *
 *    Copyright (c) 2014-2020,2022-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/binaryoutput.h"

#include <cstdint>
#include <filesystem>
#include <string>

#include "smash/action.h"
#include "smash/clock.h"
#include "smash/config.h"

namespace smash {

static auto get_list_of_binary_quantities(const std::string &content,
                                          const std::string &format,
                                          const OutputParameters &parameters);

static auto get_binary_filename(const std::string &content,
                                const std::vector<std::string> &quantities) {
  std::string filename = content;
  if (content == "Particles" || content == "Collisions") {
    std::transform(filename.begin(), filename.end(), filename.begin(),
                   [](unsigned char c) { return std::tolower(c); });
  } else if (content == "Photons" || content == "Dileptons") {
    // Nothing to be done here
  } else if (content == "Initial_Conditions") {
    filename = "SMASH_IC";
  } else {
    throw std::invalid_argument(
        "Unknown content to get the binary output filename.");
  }
  if (quantities == OutputDefaultQuantities::oscar2013) {
    filename += "_oscar2013";
  } else if (quantities == OutputDefaultQuantities::oscar2013extended) {
    filename += "_oscar2013_extended";
  } else {
    filename += "_custom";
  }
  return filename + ".bin";
}

/*!\Userguide
 * \page doxypage_output_binary
 * SMASH supports a binary output version similar to the OSCAR 2013 standard.
 * It is faster to read and write and theoretically needs less disk space.
 * However, currently in ASCII OSCAR 2013 only 5 digits after the comma are
 * written for any real number, while the binary saves the whole double
 * (16 digits). By accident, this makes the sizes of the binary output files
 * approximately the same as the OSCAR ASCII files.
 * **The binary format follows the general block structure of the OSCAR
 * format:** \ref doxypage_output_oscar. However, for the binary format,
 * the data type specification is stricter. The types used for the output are 4
 * bytes signed integers, 8 bytes doubles and 1 byte chars.
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
 * \li Format version is an integer number, currently it is 10.
 * \li Format variant is an integer number:
 *     \li 0 for quantities corresponding to OSCAR 2013 format;
 *     \li 1 for quantities corresponding to OSCAR 2013 extended format;
 *     \li 2 for custom list of quantities.
 * \li len is the length of smash version string
 * \li smash_version is len chars that give information about the SMASH version.
 *
 * **Output block header**\n
 * At start of event, end of event or any other particle output:
 * \code
 * char    int32_t       int32_t     uint32_t
 * 'p'  event_number ensemble_number n_part_lines
 * \endcode
 * \li \key event_number: Number of the event, starting with 0.
 * \li \key ensemble_number: Number of the ensemble, starting with 0.
 * \li \c n_part_lines is the number of particle lines in the block that follows
 *
 * At interaction:
 * \code
 * char uint32_t uint32_t double  double  double  uint32_t
 * 'i'  nin      nout     density xsection partial_xsection process_type
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
 *     double        int32_t         int32_t        int32_t         int32_t
 * </div>
 * <div class="line">
 * t x y z mass p0 px py pz pdg ID charge ncoll form_time xsecfac
 * proc_id_origin proc_type_origin time_last_coll pdg_mother1 pdg_mother2
 * baryon_number strangeness
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
 * \li \key baryon_number: Baryon number of the particle. 1 for baryons, -1 for
 * anti-baryons and 0 for mesons.
 *
 * **Custom Particle line**
 *
 * Similar to the \ref doxypage_output_ascii "ASCII format", the binary format
 * also supports custom quantities for particle lines. An example of particle
 * quantities is shown below:
 * \verbatim
     Output:
       Particles:
           Format:     ["Binary"]
           Quantities: ["p0", "pz", "pdg", "charge"]
   \endverbatim
 * Here, the particle data will be serialized in the same order as they appear
 * in the Quantities list. The \ref doxypage_output_ascii "ASCII format table"
 * contains the types of the quantities written in the file and be able to
 * correctly read the output e.g. in an analysis software.
 *
 * \attention If a custom binary format is used, there is no way to know which
 * quantities were stored from the output file. It is the user's responsibility
 * to keep track of this information in their projects.
 *
 * **Event end line**
 *
 * \code
 * char   int32_t       int32_t          double      char
 * 'f' event_number ensemble_number impact_parameter empty
 * \endcode
 * Where
 * \li \key event_number: Number of the event, starting with 0.
 * \li \key ensemble_number: Number of the ensemble, starting with 0.
 * \li \key impact_parameter: Impact parameter [fm] of the collision in case of
 * a collider setup, 0.0 otherwise.
 * \li \key empty: 0 if there was an interaction between the projectile
 * and the target, 1 otherwise. For non-collider setups, this is always 0.
 *
 * <h2> %Particles output </h2>
 *
 * The name of particles output file depends on its content:
 *  \li \c particles_custom.bin &rarr; this is the default;
 *  \li \c particles_oscar2013.bin &rarr;
 *      if the list of quantities corresponds to the OSCAR2013 format;
 *  \li \c particles_oscar2013_extended.bin &rarr;
 *      if the list of quantities corresponds to the extended OSCAR2013 format.
 *
 * The output file contains the current particle list at specific moments of
 * time. Every moment of time is written as a \c 'p' block. For options of this
 * output see the corresponding \ref input_output_content_specific_
 * "content-specific output options".
 *
 * <h2> Collisions output </h2>
 *
 * The name of collisions output file depends on its content:
 *  \li \c collisions_custom.bin &rarr; this is the default;
 *  \li \c collisions_oscar2013.bin &rarr;
 *      if the list of quantities corresponds to the OSCAR2013 format;
 *  \li \c collisions_oscar2013_extended.bin &rarr;
 *      if the list of quantities corresponds to the extended OSCAR2013 format.
 *
 * It contains interactions (collisions, decays, box wall crossings) and
 * optionally the initial and final configuration. The interactions are written
 * in computational frame time-ordered fashion, in \c 'i' blocks, which contains
 * the information of the incoming and the outgoing particles of each reaction
 * written in the 'incoming' and 'outgoing' blocks respectively.
 * Initial and final states are written as \c 'p' blocks. The process IDs
 * indicating the types of the reaction, such as resonance decay,
 * elastic scattering, soft string process, hard string process, etc.,
 * are written in the 'process_type' blocks. For options of this output see the
 * \ref input_output_content_specific_ "content-specific output options".
 *
 * See also \ref doxypage_output_collisions_box_modus.
 **/

BinaryOutputBase::BinaryOutputBase(const std::filesystem::path &path,
                                   const std::string &mode,
                                   const std::string &name,
                                   const std::vector<std::string> &quantities)
    : OutputInterface(name), file_{path, mode}, formatter_(quantities) {
  if (quantities.empty()) {
    throw std::invalid_argument(
        "Empty quantities list passed to 'BinaryOutputBase' constructor.");
  }
  std::fwrite("SMSH", 4, 1, file_.get());  // magic number
  write(format_version_);                  // file format version number
  std::uint16_t format_variant{};
  if (quantities == OutputDefaultQuantities::oscar2013) {
    format_variant = 0;
  } else if (quantities == OutputDefaultQuantities::oscar2013extended) {
    format_variant = 1;
  } else {
    format_variant = format_custom_;
  }
  write(format_variant);
  write(SMASH_VERSION);
}

// write functions:
void BinaryOutputBase::write(const char c) {
  std::fwrite(&c, sizeof(char), 1, file_.get());
}
void BinaryOutputBase::write(const ToBinary::type &chunk) {
  std::fwrite(chunk.data(), sizeof(char), chunk.size(), file_.get());
}

void BinaryOutputBase::write(const std::string &s) {
  const auto size = smash::numeric_cast<uint32_t>(s.size());
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
  write_in_chunk<ToBinary>(
      particles, formatter_,
      [this](const ToBinary::type &buf) { this->write(buf); });
}

void BinaryOutputBase::write(const ParticleList &particles) {
  write_in_chunk<ToBinary>(
      particles, formatter_,
      [this](const ToBinary::type &buf) { this->write(buf); });
}

void BinaryOutputBase::write_particledata(const ParticleData &p) {
  write(formatter_.particle_line(p));
}

BinaryOutputCollisions::BinaryOutputCollisions(
    const std::filesystem::path &path, std::string name,
    const OutputParameters &out_par, const std::vector<std::string> &quantities)
    : BinaryOutputBase(path / get_binary_filename(name, quantities), "wb", name,
                       quantities),
      print_start_end_(out_par.coll_printstartend) {}

void BinaryOutputCollisions::at_eventstart(const Particles &particles,
                                           const EventLabel &event_label,
                                           const EventInfo &) {
  const char pchar = 'p';
  if (print_start_end_) {
    std::fwrite(&pchar, sizeof(char), 1, file_.get());
    write(event_label.event_number);
    write(event_label.ensemble_number);
    write(particles.size());
    write(particles);
  }
}

void BinaryOutputCollisions::at_eventend(const Particles &particles,
                                         const EventLabel &event_label,
                                         const EventInfo &event) {
  const char pchar = 'p';
  if (print_start_end_) {
    std::fwrite(&pchar, sizeof(char), 1, file_.get());
    write(event_label.event_number);
    write(event_label.ensemble_number);
    write(particles.size());
    write(particles);
  }

  // Event end line
  const char fchar = 'f';
  std::fwrite(&fchar, sizeof(char), 1, file_.get());
  write(event_label.event_number);
  write(event_label.ensemble_number);
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

BinaryOutputParticles::BinaryOutputParticles(
    const std::filesystem::path &path, std::string name,
    const OutputParameters &out_par, const std::vector<std::string> &quantities)
    : BinaryOutputBase(path / get_binary_filename(name, quantities), "wb", name,
                       quantities),
      only_final_(out_par.part_only_final) {}

void BinaryOutputParticles::at_eventstart(const Particles &particles,
                                          const EventLabel &event_label,
                                          const EventInfo &) {
  const char pchar = 'p';
  if (only_final_ == OutputOnlyFinal::No) {
    std::fwrite(&pchar, sizeof(char), 1, file_.get());
    write(event_label.event_number);
    write(event_label.ensemble_number);
    write(particles.size());
    write(particles);
  }
}

void BinaryOutputParticles::at_eventend(const Particles &particles,
                                        const EventLabel &event_label,
                                        const EventInfo &event) {
  const char pchar = 'p';
  if (!(event.empty_event && only_final_ == OutputOnlyFinal::IfNotEmpty)) {
    std::fwrite(&pchar, sizeof(char), 1, file_.get());
    write(event_label.event_number);
    write(event_label.ensemble_number);
    write(particles.size());
    write(particles);
  }

  // Event end line
  const char fchar = 'f';
  std::fwrite(&fchar, sizeof(char), 1, file_.get());
  write(event_label.event_number);
  write(event_label.ensemble_number);
  write(event.impact_parameter);
  const char empty = event.empty_event;
  write(empty);

  // Flush to disk
  std::fflush(file_.get());
}

void BinaryOutputParticles::at_intermediate_time(const Particles &particles,
                                                 const std::unique_ptr<Clock> &,
                                                 const DensityParameters &,
                                                 const EventLabel &event_label,
                                                 const EventInfo &) {
  const char pchar = 'p';
  if (only_final_ == OutputOnlyFinal::No) {
    std::fwrite(&pchar, sizeof(char), 1, file_.get());
    write(event_label.event_number);
    write(event_label.ensemble_number);
    write(particles.size());
    write(particles);
  }
}

BinaryOutputInitialConditions::BinaryOutputInitialConditions(
    const std::filesystem::path &path, std::string name,
    const std::vector<std::string> &quantities)
    : BinaryOutputBase(path / get_binary_filename(name, quantities), "wb", name,
                       quantities) {}

void BinaryOutputInitialConditions::at_eventstart(const Particles &,
                                                  const EventLabel &,
                                                  const EventInfo &) {}

void BinaryOutputInitialConditions::at_eventend(
    [[maybe_unused]] const Particles &particles, const EventLabel &event_label,
    const EventInfo &event) {
  // Event end line
  const char fchar = 'f';
  std::fwrite(&fchar, sizeof(char), 1, file_.get());
  write(event_label.event_number);
  write(event_label.ensemble_number);
  write(event.impact_parameter);
  const char empty = event.empty_event;
  write(empty);

  // Flush to disk
  std::fflush(file_.get());
}

void BinaryOutputInitialConditions::at_interaction(const Action &action,
                                                   const double) {
  if (action.get_type() == ProcessType::Fluidization ||
      action.get_type() == ProcessType::FluidizationNoRemoval) {
    const char pchar = 'p';
    std::fwrite(&pchar, sizeof(char), 1, file_.get());
    write(action.incoming_particles().size());
    write(action.incoming_particles());
  }
}

static auto get_list_of_binary_quantities(const std::string &content,
                                          const std::string &format,
                                          const OutputParameters &parameters) {
  const bool is_extended = std::invoke([&content, &parameters]() {
    if (content == "Particles")
      return parameters.part_extended;
    else if (content == "Collisions")
      return parameters.coll_extended;
    else if (content == "Dileptons")
      return parameters.dil_extended;
    else if (content == "Photons")
      return parameters.photons_extended;
    else if (content == "Initial_Conditions")
      return parameters.ic_extended;
    else
      return false;
  });
  const auto default_quantities =
      (is_extended) ? OutputDefaultQuantities::oscar2013extended
                    : OutputDefaultQuantities::oscar2013;
  if (format == "Oscar2013_bin") {
    return default_quantities;
  } else if (format == "Binary") {
    if (content == "Particles" || content == "Collisions" ||
        content == "Dileptons" || content == "Photons" ||
        content == "Initial_Conditions") {
      auto list_of_quantities = parameters.quantities.at(content);
      if (list_of_quantities.empty()) {
        return default_quantities;
      } else {
        return list_of_quantities;
      }
    } else {
      /* Note that this function should not be called with "Binary" format for
       * output contents which do not support custom binary quantities. Hence we
       * throw here to prevent such a case.*/
      throw std::invalid_argument(
          "Unknown content to get the list of quantities for binary output.");
    }
  } else {
    throw std::invalid_argument(
        "Unknown format to get the list of quantities for binary output.");
  }
}

std::unique_ptr<OutputInterface> create_binary_output(
    const std::string &format, const std::string &content,
    const std::filesystem::path &path, const OutputParameters &out_par) {
  const auto quantities =
      get_list_of_binary_quantities(content, format, out_par);
  if (content == "Particles") {
    return std::make_unique<BinaryOutputParticles>(path, content, out_par,
                                                   quantities);
  } else if (content == "Collisions" || content == "Dileptons" ||
             content == "Photons") {
    return std::make_unique<BinaryOutputCollisions>(path, content, out_par,
                                                    quantities);
  } else if (content == "Initial_Conditions") {
    return std::make_unique<BinaryOutputInitialConditions>(path, content,
                                                           quantities);
  } else {
    throw std::invalid_argument("Binary output not available for '" + content +
                                "' content.");
  }
}

}  // namespace smash
