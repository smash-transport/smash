/*
 *
 *    Copyright (c) 2013-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/listmodus.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <list>
#include <map>
#include <sstream>
#include <utility>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "smash/algorithms.h"
#include "smash/angles.h"
#include "smash/configuration.h"
#include "smash/constants.h"
#include "smash/distributions.h"
#include "smash/experimentparameters.h"
#include "smash/fourvector.h"
#include "smash/inputfunctions.h"
#include "smash/logging.h"
#include "smash/macros.h"
#include "smash/particles.h"
#include "smash/random.h"
#include "smash/threevector.h"

namespace smash {

/*!\Userguide
 * \page input_modi_list_ List
 *
 * The list modus provides a modus for hydro afterburner calculations. It takes
 * files with a list of particles in \ref oscar2013_format "Oscar 2013 format"
 * as an input. These particles are treated as a starting setup. Multiple events
 * per file are supported. The input
 * parameters are:
 *
 * \key File_Directory (string, required):\n
 * Directory for the external particle lists.
 *
 * \key File_Prefix    (string, required):\n
 * Prefix for the external particle lists file.
 *
 * \key Shift_Id (int, required):\n
 * Starting id for file_id_, i.e. the first file which is read.
 *
 * \n
 * Example: Configuring an Afterburner Simulation
 * --------------
 * The following example sets up an afterburner simulation for a set of particle
 * files located in "particle_lists_in". The files are named as
 * "event{event_id}". SMASH is run once for each event in the folder.
 * \verbatim
 Modi:
     List:
         File_Directory: "particle_lists_in"
         File_Prefix: "event"

 \endverbatim
 *
 * It might for some reason be necessary to not run SMASH starting with the
 * first file. In this case, the file_id can be shifted.
 *\verbatim
 Modi:
     List:
         Shift_Id: 10
 \endverbatim
 *
 * \n
 * Example: Structure of Input Particle File
 * --------------
 * The following example shows how an input file should be formatted:
 * <div class="fragment">
 * <div class="line"><span class="preprocessor">#!OSCAR2013 particle_lists
 * t x y z mass p0 px py pz pdg ID charge</span></div>
 * <div class="line"><span class="preprocessor">\# Units: fm fm fm fm
 * GeV GeV GeV GeV GeV none none none</span></div>
 * <div class="line"><span class="preprocessor">0.1 6.42036 1.66473 9.38499
 * 0.138 0.232871 0.116953 -0.115553 0.090303 111 0 0</span></div>
 * <div class="line"><span class="preprocessor">\# event 0 end</span></div>
 * <div class="line"><span class="preprocessor">\# event 1</span></div>
 * <div class="line"><span class="preprocessor">0.1 6.42036 1.66473 9.38499
 * 0.138 0.232871 0.116953 -0.115553 0.090303 111 0 0</span></div>
 * <div class="line"><span class="preprocessor">\# event 1 end</span></div>
 * </div>
 * It means that one \f$ \pi^0 \f$ with spatial coordinates
 * (t, x, y, z) = (0.1, 6.42036, 1.66473, 9.38499) fm and
 * and 4-momenta (p0, px, py, pz) =
 * (0.232871, 0.116953, -0.115553, 0.090303) GeV,
 * with mass = 0.138 GeV, pdg = 111, id = 0 and charge 0 will be initialized for
 * the first event (and also for the second event).
 *
 * \n
 * \note
 * SMASH is shipped with an example configuration file to set up an afterburner
 * simulation by means of the list modus. This also requires a particle list to
 * be read in. Both, the configuration file and the particle list, are located
 * in /input/list. To run SMASH with the provided example configuration and
 * particle list, execute \n
 * \n
 * \verbatim
    ./smash -i INPUT_DIR/list/config.yaml
 \endverbatim
 * \n
 * Where 'INPUT_DIR' needs to be replaced by the path to the input directory
 * ('../input', if the build directory is located in the smash
 * folder).
 */

ListModus::ListModus(Configuration modus_config, const ExperimentParameters &)
    : shift_id_(modus_config.take({"List", "Shift_Id"})) {
  std::string fd = modus_config.take({"List", "File_Directory"});
  particle_list_file_directory_ = fd;

  std::string fp = modus_config.take({"List", "File_Prefix"});
  particle_list_file_prefix_ = fp;

  event_id_ = 0;
  file_id_ = shift_id_;
}

/* console output on startup of List specific parameters */
std::ostream &operator<<(std::ostream &out, const ListModus &m) {
  out << "-- List Modus\nInput directory for external particle lists:\n"
      << m.particle_list_file_directory_ << "\n";
  return out;
}

void ListModus::backpropagate_to_same_time(Particles &particles) {
  /* (1) If particles are already at the same time - don't touch them
         AND start at the start_time_ from the config. */
  double earliest_formation_time = DBL_MAX;
  double formation_time_difference = 0.0;
  double reference_formation_time = 0.0;  // avoid compiler warning
  bool first_particle = true;
  for (const auto &particle : particles) {
    const double t = particle.position().x0();
    if (t < earliest_formation_time) {
      earliest_formation_time = t;
    }
    if (first_particle) {
      reference_formation_time = t;
      first_particle = false;
    } else {
      formation_time_difference += std::abs(t - reference_formation_time);
    }
  }
  /* (2) If particles are NOT at the same time -> anti-stream them to
         the earliest time (Note: not to the start_time_ set by config) */
  bool anti_streaming_needed = (formation_time_difference > really_small);
  start_time_ = earliest_formation_time;
  if (anti_streaming_needed) {
    for (auto &particle : particles) {
      /* for hydro output where formation time is different */
      const double t = particle.position().x0();
      const double delta_t = t - start_time_;
      const ThreeVector r =
          particle.position().threevec() - delta_t * particle.velocity();
      particle.set_4position(FourVector(start_time_, r));
      particle.set_formation_time(t);
      particle.set_cross_section_scaling_factor(0.0);
    }
  }
}

void ListModus::try_create_particle(Particles &particles, PdgCode pdgcode,
                                    double t, double x, double y, double z,
                                    double mass, double E, double px, double py,
                                    double pz) {
  constexpr int max_warns_precision = 10, max_warn_mass_consistency = 10;
  const auto &log = logger<LogArea::List>();
  try {
    ParticleData &particle = particles.create(pdgcode);
    // SMASH mass versus input mass consistency check
    if (particle.type().is_stable() &&
        std::abs(mass - particle.pole_mass()) > really_small) {
      if (n_warns_precision_ < max_warns_precision) {
        log.warn() << "Provided mass of " << particle.type().name() << " = "
                   << mass << " [GeV] is inconsistent with SMASH value = "
                   << particle.pole_mass() << ". Forcing E = sqrt(p^2 + m^2)"
                   << ", where m is SMASH mass.";
        n_warns_precision_++;
      } else if (n_warns_precision_ == max_warns_precision) {
        log.warn(
            "Further warnings about SMASH mass versus input mass"
            " inconsistencies will be suppressed.");
        n_warns_precision_++;
      }
      particle.set_4momentum(mass, ThreeVector(px, py, pz));
    }
    particle.set_4momentum(FourVector(E, px, py, pz));
    // On-shell condition consistency check
    if (std::abs(particle.momentum().sqr() - mass * mass) > really_small) {
      if (n_warns_mass_consistency_ < max_warn_mass_consistency) {
        log.warn() << "Provided 4-momentum " << particle.momentum() << " and "
                   << " mass " << mass << " do not satisfy E^2 - p^2 = m^2."
                   << " This may originate from the lack of numerical"
                   << " precision in the input. Setting E to sqrt(p^2 + m^2).";
        n_warns_mass_consistency_++;
      } else if (n_warns_mass_consistency_ == max_warn_mass_consistency) {
        log.warn(
            "Further warnings about E != sqrt(p^2 + m^2) will"
            " be suppressed.");
        n_warns_mass_consistency_++;
      }
      particle.set_4momentum(mass, ThreeVector(px, py, pz));
    }
    // Set spatial coordinates, they will later be backpropagated if needed
    particle.set_4position(FourVector(t, x, y, z));
    particle.set_formation_time(t);
    particle.set_cross_section_scaling_factor(1.0);
  } catch (ParticleType::PdgNotFoundFailure) {
    log.warn() << "SMASH does not recognize pdg code " << pdgcode
               << " loaded from file. This particle will be ignored.\n";
  }
}

/* initial_conditions - sets particle data for @particles */
double ListModus::initial_conditions(Particles *particles,
                                     const ExperimentParameters &) {
  const auto &log = logger<LogArea::List>();
  std::string particle_list = next_event_();

  for (const Line &line : line_parser(particle_list)) {
    std::istringstream lineinput(line.text);
    double t, x, y, z, mass, E, px, py, pz;
    int id, charge;
    std::string pdg_string;
    lineinput >> t >> x >> y >> z >> mass >> E >> px >> py >> pz >>
        pdg_string >> id >> charge;
    if (lineinput.fail()) {
      throw LoadFailure(
          build_error_string("While loading external particle lists data:\n"
                             "Failed to convert the input string to the "
                             "expected data types.",
                             line));
    }
    PdgCode pdgcode(pdg_string);
    log.debug("Particle ", pdgcode, " (x,y,z)= (", x, ", ", y, ", ", z, ")");

    // Charge consistency check
    if (pdgcode.charge() != charge) {
      log.error() << "Charge of pdg = " << pdgcode << " != " << charge;
      throw std::invalid_argument("Inconsistent input (charge).");
    }
    try_create_particle(*particles, pdgcode, t, x, y, z, mass, E, px, py, pz);
  }
  backpropagate_to_same_time(*particles);
  event_id_++;

  return start_time_;
}

bf::path ListModus::file_path_(const int file_id) {
  const auto &log = logger<LogArea::List>();
  std::stringstream fname;
  fname << particle_list_file_prefix_ << file_id;

  const bf::path default_path = bf::absolute(particle_list_file_directory_);

  const bf::path fpath = default_path / fname.str();

  log.debug() << fpath.filename().native() << '\n';

  if (!bf::exists(fpath)) {
    log.fatal() << fpath.filename().native() << " does not exist! \n"
                << "\n Usage of smash with external particle lists:\n"
                << "1. Put the external particle lists in file \n"
                << "File_Directory/File_Prefix{id} where {id} "
                << "traversal [Shift_Id, Nevent-1]\n"
                << "2. Particles info: t x y z mass p0 px py pz"
                << " pdg ID charge\n"
                << "in units of: fm fm fm fm GeV GeV GeV GeV GeV"
                << " none none none\n";
    throw std::runtime_error("External particle list does not exist!");
  }

  return fpath;
}

std::string ListModus::next_event_() {
  const auto &log = logger<LogArea::List>();

  const bf::path fpath = file_path_(file_id_);
  bf::ifstream ifs{fpath};
  ifs.seekg(last_read_position_);

  if (!file_has_events_(fpath, last_read_position_)) {
    // current file out of events. get next file and call this function
    // recursively.
    file_id_++;
    last_read_position_ = 0;
    ifs.close();
    return next_event_();
  }

  // read one event. events marked by line # event end i in case of Oscar
  // output. Assume one event per file for all other output formats
  std::string event_string;
  const std::string needle = "end";
  std::string line;
  while (getline(ifs, line)) {
    if (line.find(needle) == std::string::npos) {
      event_string += line + "\n";
    } else {
      break;
    }
  }

  if (!ifs.eof() && (ifs.fail() || ifs.bad())) {
    log.fatal() << "Error while reading " << fpath.filename().native();
    throw std::runtime_error("Error while reading external particle list");
  }
  // save position for next event read
  last_read_position_ = ifs.tellg();
  ifs.close();

  return event_string;
}

bool ListModus::file_has_events_(bf::path filepath,
                                 std::streampos last_position) {
  const auto &log = logger<LogArea::List>();
  bf::ifstream ifs{filepath};
  std::string line;

  // last event read read at end of file. we know this because errors are
  // handled in next_event
  if (last_position == -1) {
    return false;
  }
  ifs.seekg(last_position);
  // skip over comment lines, assume that a max. of four consecutive comment
  // lines can occur
  int skipped_lines = 0;
  const int max_comment_lines = 4;
  while (std::getline(ifs, line) && line[0] != '#' &&
         skipped_lines++ < max_comment_lines) {
  }

  if (ifs.eof()) {
    return false;
  }

  if (!ifs.good()) {
    log.fatal() << "Error while reading " << filepath.filename().native();
    throw std::runtime_error("Error while reading external particle list");
  }

  ifs.close();
  return true;
}

}  // namespace smash
