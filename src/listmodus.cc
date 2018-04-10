/*
 *
 *    Copyright (c) 2013-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/listmodus.h"

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

#include "include/algorithms.h"
#include "include/angles.h"
#include "include/configuration.h"
#include "include/constants.h"
#include "include/distributions.h"
#include "include/experimentparameters.h"
#include "include/fourvector.h"
#include "include/inputfunctions.h"
#include "include/logging.h"
#include "include/macros.h"
#include "include/particles.h"
#include "include/random.h"
#include "include/threevector.h"

namespace smash {

/*!\Userguide
 * \page input_modi_list_ List
 *
 * <b> Purpose: Provides a modus for hydro afterburner calculations </b>
 *
 * Takes files with list of particles in
 * \ref oscar2013_format "Oscar 2013 format" as an input. These
 * particles are treated as a starting setup.
 *
 * To select the List modus:
 * \code
 * General:
 *     Modus:  List
 * \endcode
 *
 * ### Input options:
 *
 * \key File_Directory (string, required):\n
 * Directory for the external particle lists
 *
 * \key File_Prefix    (string, required):\n
 * Prefix for the external particle lists file
 *
 * \key Start_Time (double, required):\n
 * Starting time of List calculation.
 *
 * \key Shift_Id (int, required):\n
 * Starting id for event_id_
 *
 * ### Example of configuration:
 *
 * \code
 *  List:
 *      # path of external particle list == File_Directory/File_Prefix{eventid}
 *      File_Directory: "particle_lists_in"
 *      File_Prefix: "event"
 *
 *      # starting number of file_id in event-by-event simulation
 *      Shift_Id: 0
 *
 *      # start time of simulation
 *      Start_Time: 0.0
 * \endcode
 *
 * ### Example of input in particle_lists_in/event{id}
 * <div class="fragment">
 * <div class="line"><span class="preprocessor">#!OSCAR2013 particle_lists
 * t x y z mass p0 px py pz pdg ID charge</span></div>
 * <div class="line"><span class="preprocessor">\# Units: fm fm fm fm
 * GeV GeV GeV GeV GeV none none none</span></div>
 * <div class="line"><span class="preprocessor">0.1 6.42036 1.66473 9.38499
 * 0.138 0.232871 0.116953 -0.115553 0.090303 111 0 0</span></div>
 * </div>
 * It means that one \f$ \pi^0 \f$ with spatial coordinates\n
 * (t, x, y, z) = (0.1, 6.42036, 1.66473, 9.38499) fm and\n
 * and 4-momenta (p0, px, py, pz) =
 * (0.232871, 0.116953, -0.115553, 0.090303) GeV,\n
 * with mass = 0.138 GeV, pdg = 111, id = 0 and charge 0 will be initialized.
 *
 * ### WARNING:
 * Currently only one event per file is supported. Having more than one
 * event per file will lead to undefined behavior.
 */

ListModus::ListModus(Configuration modus_config, const ExperimentParameters &)
    : start_time_(modus_config.take({"List", "Start_Time"})),
      shift_id_(modus_config.take({"List", "Shift_Id"})) {
  std::string fd = modus_config.take({"List", "File_Directory"});
  particle_list_file_directory_ = fd;

  std::string fp = modus_config.take({"List", "File_Prefix"});
  particle_list_file_prefix_ = fp;

  event_id_ = 0;
  file_id_ = shift_id_;
}

/* console output on startup of List specific parameters */
std::ostream &operator<<(std::ostream &out, const ListModus &m) {
  out << "\nStarting time for List calculation: " << m.start_time_ << '\n';
  out << "\nInput directory for external particle lists:"
      << m.particle_list_file_directory_ << "\n";
  return out;
}

std::pair<bool, double> ListModus::check_formation_time_(
    const std::string &particle_list) {
  double earliest_formation_time = DBL_MAX;
  double formation_time_difference = 0.0;
  double reference_formation_time = 0.0;  // avoid compiler warning
  for (const Line &line : line_parser(particle_list)) {
    std::istringstream lineinput(line.text);
    double t;
    lineinput >> t;
    if (t < earliest_formation_time) {
      earliest_formation_time = t;
    }

    if (line.number == 0) {
      reference_formation_time = t;
    } else {
      formation_time_difference += std::abs(t - reference_formation_time);
    }
  }

  bool anti_streaming_needed =
      (formation_time_difference > really_small) ? true : false;
  return std::make_pair(anti_streaming_needed, earliest_formation_time);
}

/* initial_conditions - sets particle data for @particles */
double ListModus::initial_conditions(Particles *particles,
                                     const ExperimentParameters &) {
  const auto &log = logger<LogArea::List>();

  std::string particle_list = next_event_();

  auto check = check_formation_time_(particle_list);
  bool anti_streaming_needed = std::get<0>(check);
  start_time_ = std::get<1>(check);

  constexpr int max_warns_precision = 10, max_warn_mass_consistency = 10;

  for (const Line &line : line_parser(particle_lists)) {
    std::istringstream lineinput(line.text);
    double t, x, y, z, mass, E, px, py, pz;
    int id, charge;
    PdgCode pdgcode;
    lineinput >> t >> x >> y >> z >> mass >> E >> px >> py >> pz >> pdgcode >>
        id >> charge;

    if (lineinput.fail()) {
      throw LoadFailure(
          build_error_string("While loading external particle lists data:\n"
                             "Failed to convert the input string to the "
                             "expected data types.",
                             line));
    }

    log.debug("Particle ", pdgcode, " (x,y,z)= (", x, ", ", y, ", ", z, ")");

    try {
      ParticleData &particle = particles->create(pdgcode);
      // Charge consistency check
      if (pdgcode.charge() != charge) {
        log.error() << "Charge of pdg = " << pdgcode << " != " << charge;
        throw std::invalid_argument("Inconsistent input (charge).");
      }
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
          log.warn("Further warnings about SMASH mass versus input mass"
                   " inconsistencies will be suppressed.");
         n_warns_precision_++;
        }
        particle.set_4momentum(mass, ThreeVector(px, py, pz));
      }
      particle.set_4momentum(FourVector(E, px, py, pz));
      // On-shell condition consistency check
      if (std::abs(particle.momentum().sqr() - mass*mass) > really_small) {
        if (n_warns_mass_consistency_ < max_warn_mass_consistency) {
          log.warn() << "Provided 4-momentum " << particle.momentum() << " and "
                    << " mass " << mass << " do not satisfy E^2 - p^2 = m^2."
                    << " This may originate from the lack of numerical"
                    << " precision in the input. Setting E to sqrt(p^2 + m^2).";
          n_warns_mass_consistency_++;
        } else if (n_warns_mass_consistency_ == max_warn_mass_consistency) {
          log.warn("Further warnings about E != sqrt(p^2 + m^2) will"
                   " be suppressed.");
         n_warns_mass_consistency_++;
        }
        particle.set_4momentum(mass, ThreeVector(px, py, pz));
      }
      if (anti_streaming_needed) {
        /* for hydro output where formation time is different */
        double delta_t = t - start_time_;
        FourVector start_timespace =
            FourVector(t, x, y, z) - delta_t * FourVector(E, px, py, pz) / E;
        particle.set_4position(start_timespace);
        particle.set_formation_time(t);
        particle.set_cross_section_scaling_factor(0.0);
      } else {
        /* for smash output where formation time is the same */
        particle.set_4position(FourVector(t, x, y, z));
        particle.set_formation_time(t);
        particle.set_cross_section_scaling_factor(1.0);
      }
    } catch (ParticleType::PdgNotFoundFailure) {
      log.warn() << "While loading external particle lists data, "
                 << "PDG code not found for the particle:\n"
                 << line.text << std::endl;
    }
  }

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

  constexpr char comment_token = '#';

  const bf::path fpath = file_path_(file_id_);
  bf::ifstream ifs{fpath};
  ifs.seekg(last_read_position_);

  if (last_read_position_ == 0) {
    // skip header (4 lines starting with #)
    std::string tmp;
    while (std::getline(ifs, tmp) && ifs.peek() == comment_token)
      ;
  }

  if (!file_has_events_(fpath, last_read_position_)) {
    // current file out of events. get next file and call this function recursively.
    file_id_++;
    last_read_position_ = 0;
    ifs.close();
    return next_event_();
  }

  // read one event. events marked by line # event end i 
  std::string event_string;
  const std::string needle = "end";
  std::string line;
  while (getline(ifs, line)) {
    if (line.find(needle) == std::string::npos) {
      event_string += line + "\n";
    } else
      break;
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

  ifs.seekg(last_position);
  // skip over comment lines
  while (std::getline(ifs, line) && line[0] != '#')
    ;
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
